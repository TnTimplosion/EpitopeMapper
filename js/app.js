/* Wiring: files in, analysis, 3D view, table, exports.
 *
 * One pipeline drives everything. Any change to a setting, an exclusion or a file re-enters
 * recompute(), which runs the chunked SASA job, rebuilds the per-residue rows and the patch sums,
 * and then restyles the view and refills the table from the same rows.
 */
(function (global) {
  'use strict';

  const Chem = global.Chem, PDB = global.PDB, SASA = global.SASA, Colour = global.Colour,
        Viewer = global.Viewer, Table = global.Table;
  const $ = (id) => document.getElementById(id);

  const S = {
    st: null,                  // parsed structure
    fileName: '',
    membrane: null,            // active bilayer definition, or null
    detected: null,            // bilayer found in the loaded files
    ppmName: '',
    excluded: null,            // Uint8Array over residues
    picked: new Set(),
    focus: null,
    hovered: null,
    rows: [],
    byRid: [],
    job: null,
    elec: null,                // { phi: per-atom potential, sites: the charges behind it }
    token: 0,
    colourVersion: 0,
    chainFilter: null,         // Set of chain ids, or null for all
    ctx: {},
    scales: {},
    lastMs: 0,
    viewerReady: false,
  };

  let timer = null;

  // --------------------------------------------------------------------------------- utilities

  function fmt(v, d) {
    return (v === null || v === undefined || !Number.isFinite(v)) ? '&mdash;' : v.toFixed(d);
  }

  function fmtCharge(q) {
    if (!Number.isFinite(q)) return '&mdash;';
    if (Math.abs(q) < 1e-9) return '0';
    const a = Math.abs(q);
    return (q > 0 ? '+' : '&minus;') + (Number.isInteger(a) ? a.toFixed(0) : a.toFixed(1));
  }

  function message(el, html, kind) {
    if (!html) { el.hidden = true; el.innerHTML = ''; return; }
    el.hidden = false;
    el.className = 'msg' + (kind ? ' ' + kind : '');
    el.innerHTML = html;
  }

  function busy(show, text, frac) {
    const box = $('busy');
    box.hidden = !show;
    if (text) $('busyText').innerHTML = text;
    $('busyBar').style.width = Math.round(100 * Math.min(1, Math.max(0, frac || 0))) + '%';
  }

  function download(name, text, mime) {
    const blob = new Blob([text], { type: mime || 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  function percentile(values, p) {
    if (!values.length) return NaN;
    const v = values.slice().sort((a, b) => a - b);
    const i = Math.min(v.length - 1, Math.max(0, Math.round((v.length - 1) * p)));
    return v[i];
  }

  // ----------------------------------------------------------------------------- settings read

  function analysisOpts() {
    return {
      probe: parseFloat($('probe').value) || SASA.PROBE,
      nPoints: parseInt($('points').value, 10) || 100,
      includeH: $('incH').checked,
      includeHet: $('incHet').checked,
      surfaceRSA: parseFloat($('rsaCut').value),
      patchR: parseFloat($('patchR').value) || 10,
      scale: $('scaleSel').value || 'kd',
      headgroup: parseFloat($('memHead').value) || 0,
    };
  }

  function candOpts() {
    return {
      proteinOnly: $('aaOnly').checked,
      skipMembrane: !!S.membrane && $('memSkip').checked,
      keepHeadgroup: $('memKeepHead').checked,
    };
  }

  /** What to call each face of the bilayer, once the user has said which way round it is. */
  function sideNames() {
    const side = $('memSide').value;
    return {
      plus: side === 'out' ? 'extracellular' : side === 'in' ? 'cytoplasmic' : '+z side',
      minus: side === 'out' ? 'cytoplasmic' : side === 'in' ? 'extracellular' : '−z side',
      labelled: side !== 'unknown',
    };
  }

  function zoneLabels() {
    const n = sideNames();
    const plus = n.labelled ? n.plus : 'aqueous, +z';
    const minus = n.labelled ? n.minus : 'aqueous, &minus;z';
    return {
      'core': 'bilayer', 'head+': 'headgroup +z', 'head-': 'headgroup &minus;z',
      'side+': plus, 'side-': minus, plusLabel: plus, minusLabel: minus,
    };
  }

  // --------------------------------------------------------------------------- membrane set-up

  function structureExtent() {
    const a = S.st.atoms;
    let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
    for (const at of a) {
      if (at.x < xMin) xMin = at.x; if (at.x > xMax) xMax = at.x;
      if (at.y < yMin) yMin = at.y; if (at.y > yMax) yMax = at.y;
    }
    return { xMin: xMin - 4, xMax: xMax + 4, yMin: yMin - 4, yMax: yMax + 4 };
  }

  function updateMembrane() {
    const mode = $('memMode').value;
    const manual = mode === 'manual';
    $('memThick').disabled = !manual;
    $('memCentre').disabled = !manual;

    if (mode === 'off' || !S.st) {
      S.membrane = null;
    } else if (manual) {
      const t = parseFloat($('memThick').value) || 30;
      const c = parseFloat($('memCentre').value) || 0;
      S.membrane = {
        zLo: c - t / 2, zHi: c + t / 2, thickness: t, centre: c, normal: 'z',
        extent: structureExtent(), source: 'set by hand',
      };
    } else {
      S.membrane = S.detected ? Object.assign({}, S.detected) : null;
      if (S.membrane && !S.membrane.extent) S.membrane.extent = structureExtent();
    }

    const st = $('memStatus');
    if (!S.membrane) {
      if (S.detected && !S.st) {
        st.innerHTML = 'Bilayer recorded from the dummy atoms: ' + S.detected.thickness.toFixed(1)
          + '&nbsp;&Aring; thick, boundaries at z = ' + S.detected.zLo.toFixed(1) + ' and '
          + S.detected.zHi.toFixed(1) + '&nbsp;&Aring;. Load the structure to use it.';
      } else {
        st.innerHTML = S.detected
          ? 'Bilayer ignored. The whole surface is treated as solvent exposed.'
          : 'No membrane defined. Load a PPM 3.0 output, or set the bilayer by hand, if the target '
            + 'is a membrane protein.';
      }
    } else {
      const m = S.membrane;
      let html = 'Bilayer ' + m.thickness.toFixed(1) + '&nbsp;&Aring; thick, normal along z, '
        + 'boundaries at z = ' + m.zLo.toFixed(1) + ' and ' + m.zHi.toFixed(1) + '&nbsp;&Aring; ('
        + m.source + ').';
      if (S.st) {
        let inside = false;
        for (const a of S.st.atoms) {
          if (a.z >= m.zLo && a.z <= m.zHi) { inside = true; break; }
        }
        if (!inside) {
          html += ' <span class="warn">No residue lies inside the slab &mdash; check that the '
            + 'structure and the PPM output are the same, aligned coordinates.</span>';
        }
      }
      st.innerHTML = html;
    }
    const on = !!S.membrane;
    $('memOcclude').disabled = !on;
    $('memSkip').disabled = !on;
    $('memKeepHead').disabled = !on || !$('memSkip').checked;
    $('memShow').disabled = !on;
    $('memOpacity').disabled = !on;
    updateSideFilter();
  }

  /**
   * The membrane-face filter in the table toolbar. It only means anything once a bilayer is
   * defined, so it stays out of the way until then, and its options take the names the user has
   * given the two faces.
   */
  function updateSideFilter() {
    const group = $('sideGroup');
    const sel = $('sideSel');
    const on = !!S.membrane;
    group.hidden = !on;
    if (!on) {
      sel.value = 'all';              // a hidden filter must not go on hiding rows
      return;
    }
    const n = sideNames();
    const cap = (s) => s.charAt(0).toUpperCase() + s.slice(1);
    sel.options[1].textContent = cap(n.plus) + ' only';
    sel.options[2].textContent = cap(n.minus) + ' only';
  }

  // ---------------------------------------------------------------------------- file handling

  function readFile(file) {
    return new Promise((resolve, reject) => {
      const fr = new FileReader();
      fr.onload = () => resolve(String(fr.result || ''));
      fr.onerror = () => reject(new Error('"' + file.name + '" could not be read.'));
      fr.readAsText(file);
    });
  }

  const DUM_RE = /^(?:ATOM  |HETATM).{11}DUM/m;

  async function handleFile(file, forcedKind) {
    try {
      message($('fileMsg'), '');
      busy(true, 'Reading ' + file.name + '&hellip;', 0.2);
      const text = await readFile(file);
      // A PPM file may hold nothing but the two planes of dummy atoms; that is still a bilayer,
      // whichever slot it was given to.
      const parsed = PDB.parse(text, file.name, { allowEmpty: true });
      const detected = PDB.detectMembrane(parsed.dummies);
      let kind = forcedKind || (DUM_RE.test(text) ? 'ppm' : 'struct');
      if (!parsed.atoms.length) kind = 'ppm';

      if (kind === 'ppm') {
        if (!detected) {
          throw new Error('"' + file.name + '" holds no usable DUM dummy atoms. A PPM 3.0 result '
            + 'file carries two planes of them, which is what defines the bilayer.');
        }
        S.detected = detected;
        S.ppmName = file.name;
        $('ppmInfo').textContent = file.name;
        $('ppmInfo').classList.add('ok');
        $('memMode').value = 'auto';
        $('memThick').value = detected.thickness.toFixed(1);
        $('memCentre').value = detected.centre.toFixed(1);

        if (parsed.atoms.length) {
          // The PPM output is rotated into the membrane frame, so its coordinates are the ones the
          // slab refers to: use them for the analysis as well.
          await setStructure(parsed, file.name);
          message($('fileMsg'), 'Loaded the PPM 3.0 output as the structure too: its coordinates '
            + 'are rotated so the bilayer normal is z, which is the frame the dummy atoms belong to.');
        } else if (!S.st) {
          updateMembrane();
          message($('fileMsg'), 'Bilayer recorded. Load the structure to start the analysis.');
        } else {
          updateMembrane();
          if (S.membrane && S.viewerReady) Viewer.orientMembrane();
          recompute();
        }
        busy(false);
        return;
      }

      $('structInfo').textContent = file.name;
      $('structInfo').classList.add('ok');
      if (detected) {
        S.detected = detected;
        $('memThick').value = detected.thickness.toFixed(1);
        $('memCentre').value = detected.centre.toFixed(1);
      }
      await setStructure(parsed, file.name);
      busy(false);
    } catch (err) {
      busy(false);
      message($('fileMsg'), String(err.message || err), 'err');
    }
  }

  async function setStructure(parsed, name) {
    S.st = parsed;
    S.fileName = name;
    S.excluded = new Uint8Array(parsed.residues.length);
    S.picked = new Set();
    S.focus = null;
    S.hovered = null;
    S.chainFilter = null;

    $('empty').hidden = true;
    $('readout').hidden = true;

    let note = '';
    if (S.viewerReady) {
      const info = Viewer.load(parsed);
      if (info.matched === 'failed') {
        note = ' <span class="warn">The 3D view could not be matched to the parsed atoms, so '
          + 'highlighting may be wrong.</span>';
      } else if (info.matched === 'key') {
        note = ' Atoms were matched to the view by name.';
      }
    }

    const c = parsed.counts;
    $('structStats').innerHTML = c.atoms.toLocaleString() + ' atoms &middot; '
      + c.residues.toLocaleString() + ' residues &middot; '
      + (parsed.chains.length === 1 ? 'chain ' : 'chains ')
      + parsed.chains.map((x) => x.id).join(', ')
      + (c.waters ? ' &middot; ' + c.waters + ' waters dropped' : '')
      + (c.hydrogens ? ' &middot; ' + c.hydrogens + ' hydrogens' : '')
      + (c.dummies ? ' &middot; ' + c.dummies + ' dummy atoms' : '')
      + (c.hetResidues ? ' &middot; ' + c.hetResidues + ' het groups' : '')
      + note;

    if (parsed.warnings.length) {
      message($('fileMsg'), parsed.warnings.join('<br>'));
    }

    buildChainChips();
    updateMembrane();
    if (S.membrane && S.viewerReady) Viewer.orientMembrane();
    if (parsed.counts.atoms > 20000 && $('surfSel').value !== 'none') {
      $('surfSel').value = 'none';
    }
    recompute();
  }

  // -------------------------------------------------------------------------------- the engine

  function schedule(delay) {
    if (timer) clearTimeout(timer);
    timer = setTimeout(() => { timer = null; recompute(); }, delay === undefined ? 180 : delay);
  }

  function recompute() {
    if (!S.st) return;
    const o = analysisOpts();
    const token = ++S.token;
    const t0 = performance.now();
    const occlude = !!S.membrane && $('memOcclude').checked;

    S.job = SASA.start(S.st, {
      probe: o.probe,
      nPoints: o.nPoints,
      includeH: o.includeH,
      includeHet: o.includeHet,
      excluded: S.excluded,
      slab: occlude ? S.membrane : null,
    });

    busy(true, 'Rolling a ' + o.probe.toFixed(2) + '&nbsp;&Aring; probe over '
      + S.job.nActive.toLocaleString() + ' atoms&hellip;', 0);

    const pump = () => {
      if (token !== S.token) return;                 // a newer run took over
      if (!S.job.step || S.job.step(25)) {
        S.lastMs = performance.now() - t0;
        finish(o);
        busy(false);
        return;
      }
      busy(true, null, S.job.progress);
      setTimeout(pump, 0);
    };
    pump();
  }

  function finish(o) {
    const tStart = performance.now();
    S.rows = SASA.residueMetrics(S.st, S.job, {
      scale: o.scale,
      surfaceRSA: o.surfaceRSA,
      excluded: S.excluded,
      slab: S.membrane,
      headgroup: o.headgroup,
    });
    SASA.patchMetrics(S.rows, o.patchR, candOpts());

    // Electrostatics: per atom for the surface, area-weighted per residue for the table.
    S.elec = Elec.compute(S.st, S.job, { ionicStrength: Elec.DEFAULT_IONIC, probe: o.probe });
    const espRes = Elec.perResidue(S.st, S.job, S.elec.phi);
    for (const r of S.rows) r.esp = espRes[r.rid];

    S.byRid = new Array(S.st.residues.length);
    for (const r of S.rows) S.byRid[r.rid] = r;

    // Domains for the colour ramps, from the data rather than fixed limits.
    const surf = S.rows.filter((r) => r.surface && !r.excluded);
    const bs = S.rows.map((r) => r.b).filter(Number.isFinite);
    const chainIndex = {};
    S.st.chains.forEach((c, i) => { chainIndex[c.id] = i; });
    let bLo = bs.length ? Math.min.apply(null, bs) : 0;
    let bHi = bs.length ? Math.max.apply(null, bs) : 1;
    if (bHi - bLo < 1e-6) bHi = bLo + 1;

    // Robust upper ends for the ramps: the 95th percentile, so one huge arginine does not wash the
    // rest of the surface out. Falls back to a fixed value when nothing is exposed at all.
    const p95 = (values, fallback) => {
      const v = percentile(values, 0.95);
      return Number.isFinite(v) ? v : fallback;
    };
    S.ctx = {
      scale: o.scale,
      chains: S.st.chains,
      chainIndex: chainIndex,
      rsaMax: 100,
      sasaMax: Math.max(20, p95(surf.map((r) => r.sasa), 150)),
      patchSASAMax: Math.max(50, p95(surf.map((r) => r.patchSASA), 600)),
      patchChargeMax: Math.max(2, Math.ceil(p95(surf.map((r) => Math.abs(r.patchCharge)), 4))),
      lipidMax: Math.max(10, p95(S.rows.map((r) => r.lipid).filter((v) => v > 0), 120)),
      espMax: (() => {
        const abs = [];
        for (let i = 0; i < S.elec.phi.length; i++) {
          const v = S.elec.phi[i];
          if (Number.isFinite(v)) abs.push(Math.abs(v));
        }
        return Math.min(25, Math.max(1, p95(abs, 5)));
      })(),
      bRange: [bLo, bHi],
      bLabel: bHi <= 100.5 && bLo >= 0 && bHi > 40 ? 'B-factor or pLDDT' : 'B-factor',
      plusLabel: zoneLabels().plusLabel,
      minusLabel: zoneLabels().minusLabel,
    };
    S.scales = {
      hyd: Colour.scaleFor('hydro', S.ctx),
      rsa: Colour.scaleFor('rsa', S.ctx),
      charge: Colour.scaleFor('charge', S.ctx),
      esp: Colour.scaleFor('esp', S.ctx),
    };
    S.colourVersion++;

    const nCand = S.rows.reduce((n, r) => n + (SASA.isCandidate(r, candOpts()) ? 1 : 0), 0);
    let totalSASA = 0;
    for (const r of S.rows) if (!r.excluded) totalSASA += r.sasa;
    S.lastMs += performance.now() - tStart;
    $('analysisInfo').innerHTML = nCand + ' candidate residues &middot; '
      + Math.round(totalSASA).toLocaleString() + '&nbsp;&Aring;&sup2; accessible &middot; '
      + (S.lastMs < 1000 ? Math.round(S.lastMs) + '&nbsp;ms' : (S.lastMs / 1000).toFixed(1) + '&nbsp;s');

    refillTable();
    restyle();
    updateSummary();
    updateExclChips();
  }

  function refillTable() {
    Table.setRows(S.rows, {
      show: $('showSel').value,
      side: $('sideSel').value,
      search: $('searchBox').value,
      chains: S.chainFilter,
      picked: S.picked,
      candOpts: candOpts(),
      zoneLabels: zoneLabels(),
      hydScale: S.scales.hyd,
      rsaScale: S.scales.rsa,
      qScale: S.scales.charge,
      espScale: S.scales.esp,
      scale: $('scaleSel').value,
    });
    Table.setFocus(S.focus, false);
    $('tableEmpty').hidden = Table.visibleRows().length > 0;
  }

  function restyle() {
    if (!S.viewerReady || !Viewer.loaded()) return;
    const mode = $('colourSel').value;
    const scale = Colour.scaleFor(mode, S.ctx);
    const row = S.focus !== null ? S.byRid[S.focus] : null;
    Viewer.style({
      rows: S.rows,
      mode: mode,
      scale: scale,
      ctx: S.ctx,
      rep: $('repSel').value,
      surface: $('surfSel').value,
      opacity: parseFloat($('opacity').value),
      picked: S.picked,
      focus: S.focus,
      active: S.job ? S.job.active : null,
      // A potential is a property of a point, not of a residue, so it colours atom by atom.
      atomValues: (scale.perAtom && S.elec) ? S.elec.phi : null,
      atomScale: scale.perAtom ? scale : null,
      showExcluded: $('showExcl').checked,
      membrane: S.membrane,
      showMembrane: $('memShow').checked,
      membraneOpacity: parseFloat($('memOpacity').value),
      label: ($('showLabel').checked && row)
        ? row.label + ' ' + row.chainID + '  RSA ' + (Number.isFinite(row.rsa) ? row.rsa.toFixed(0) + '%' : '-')
        : null,
      colourKey: S.colourVersion + '|' + mode,
      surfaceSelKey: S.colourVersion,
    });
    updateLegend(mode, scale);
  }

  function updateLegend(mode, scale) {
    const box = $('legend');
    if (!S.st) { box.hidden = true; return; }
    box.hidden = false;
    $('lgTitle').innerHTML = scale.label || '';
    if (scale.categorical) {
      $('lgBar').hidden = true;
      $('lgLo').innerHTML = ''; $('lgHi').innerHTML = '';
      const cats = Colour.categories(mode, S.ctx);
      $('lgKeys').innerHTML = cats.map((c) =>
        '<span><i style="background:#' + c.colour.toString(16).padStart(6, '0') + '"></i>'
        + c.label + '</span>').join('');
    } else {
      $('lgBar').hidden = false;
      const stops = scale.invert ? scale.stops.slice().reverse() : scale.stops;
      $('lgBar').style.background = Colour.gradient(stops);
      const d = scale.domain;
      const lo = scale.loLabel, hi = scale.hiLabel;
      const num = (v) => (Math.abs(v) >= 100 ? v.toFixed(0) : v.toFixed(Math.abs(v) < 10 ? 1 : 0));
      $('lgLo').innerHTML = num(d[0]) + (lo ? ' <span class="dim">' + lo + '</span>' : '');
      $('lgHi').innerHTML = (hi ? '<span class="dim">' + hi + '</span> ' : '') + num(d[1]);
      $('lgKeys').innerHTML = '';
    }
  }

  // ----------------------------------------------------------------------------- readout panel

  function updateReadout(rid) {
    const box = $('readout');
    const r = (rid === null || rid === undefined) ? null : S.byRid[rid];
    if (!r) { box.hidden = true; return; }
    box.hidden = false;
    const sc = Chem.SCALES[$('scaleSel').value] || Chem.SCALES.kd;
    $('roName').innerHTML = r.resName + ' ' + r.resSeq + (r.iCode || '')
      + (r.res.ss ? ' <span class="ro-ss">' + r.res.ss + '</span>' : '');
    $('roChain').textContent = 'chain ' + r.chainID;
    $('roSasa').innerHTML = fmt(r.sasa, 1) + ' &Aring;&sup2;';
    $('roRsa').innerHTML = Number.isFinite(r.rsa) ? fmt(r.rsa, 1) + ' %' : '&mdash;';
    $('roHyd').innerHTML = fmt(r.hyd, 2) + (sc.units ? ' ' + sc.units : '');
    $('roCharge').innerHTML = fmtCharge(r.charge);
    $('roEsp').innerHTML = Number.isFinite(r.esp)
      ? (r.esp > 0 ? '+' : '') + fmt(r.esp, 1) + ' kcal/mol&middot;e' : '&mdash;';
    $('roPatch').innerHTML = r.patchN ? fmt(r.patchSASA, 0) + ' &Aring;&sup2;, ' + r.patchN + ' res'
      : '&mdash;';
    $('roZone').innerHTML = r.zone ? (zoneLabels()[r.zone] || r.zone) : '&mdash;';
    const notes = [];
    if (r.excluded) notes.push('excluded from the analysis');
    else if (!r.surface) notes.push('below the surface cut');
    else if (!SASA.isCandidate(r, candOpts())) notes.push('not an epitope candidate');
    if (r.lipid > 1) notes.push(fmt(r.lipid, 0) + ' &Aring;&sup2; buried in the bilayer');
    const note = $('roNote');
    note.hidden = !notes.length;
    note.innerHTML = notes.join(' &middot; ');
  }

  function updateSummary() {
    const rows = [];
    for (const rid of S.picked) if (S.byRid[rid]) rows.push(S.byRid[rid]);
    const sum = SASA.summarise(rows, $('scaleSel').value);
    $('pkN').textContent = sum.n;
    $('pkSasa').innerHTML = sum.n ? fmt(sum.sasa, 0) + ' &Aring;&sup2;' : '&mdash;';
    $('pkHyd').innerHTML = sum.n ? fmt(sum.hyd, 2) : '&mdash;';
    $('pkCharge').innerHTML = sum.n ? fmtCharge(sum.charge) : '&mdash;';
    $('pkSpan').innerHTML = sum.n ? fmt(sum.span, 1) + ' &Aring;' : '&mdash;';
    $('hotspots').value = hotspotString();
  }

  function hotspotString() {
    const rows = [];
    for (const rid of S.picked) if (S.byRid[rid]) rows.push(S.byRid[rid]);
    rows.sort((a, b) => a.chainID.localeCompare(b.chainID) || a.resSeq - b.resSeq
      || a.iCode.localeCompare(b.iCode));
    return rows.map((r) => r.tag).join(',');
  }

  // ------------------------------------------------------------------------ selection handling

  function focusResidue(rid, opts) {
    opts = opts || {};
    S.focus = rid;
    Table.setFocus(rid, opts.scroll !== false);
    updateReadout(rid);
    restyle();
    if (rid !== null && (opts.centre || $('centreClick').checked)) Viewer.centreOn(rid, false);
  }

  function togglePick(rid) {
    if (S.picked.has(rid)) S.picked.delete(rid); else S.picked.add(rid);
    S.focus = rid;
    refillTable();
    updateReadout(rid);
    updateSummary();
    restyle();
  }

  function toggleDrop(rid) {
    S.excluded[rid] = S.excluded[rid] ? 0 : 1;
    if (S.excluded[rid]) S.picked.delete(rid);
    recompute();
  }

  function setExcluded(rids, on) {
    for (const rid of rids) {
      S.excluded[rid] = on ? 1 : 0;
      if (on) S.picked.delete(rid);
    }
    recompute();
  }

  /** Residue selection language: "A12", "A20-28", "12-20", "B", "het". */
  function parseSelection(text) {
    const out = new Set(), bad = [];
    if (!S.st) return { rids: out, bad: bad };
    const chainIds = new Set(S.st.chains.map((c) => c.id));
    for (const tok of String(text).split(/[,;\s]+/)) {
      const t = tok.trim();
      if (!t) continue;
      if (/^het$/i.test(t)) {
        for (const r of S.st.residues) if (!r.isAA) out.add(r.rid);
        continue;
      }
      if (t.length === 1 && chainIds.has(t)) {
        for (const r of S.st.residues) if (r.chainID === t) out.add(r.rid);
        continue;
      }
      let m = /^([A-Za-z])?(-?\d+)([A-Za-z])?$/.exec(t);
      if (m) {
        const ch = m[1] || null, num = parseInt(m[2], 10), ic = m[3] || '';
        let hit = 0;
        for (const r of S.st.residues) {
          if (ch && r.chainID !== ch) continue;
          if (r.resSeq !== num) continue;
          if (ic && r.iCode.toUpperCase() !== ic.toUpperCase()) continue;
          out.add(r.rid); hit++;
        }
        if (!hit) bad.push(t);
        continue;
      }
      m = /^([A-Za-z])?(-?\d+)([A-Za-z])?-([A-Za-z])?(-?\d+)([A-Za-z])?$/.exec(t);
      if (m) {
        const ch = m[1] || m[4] || null;
        const a = parseInt(m[2], 10), b = parseInt(m[5], 10);
        const lo = Math.min(a, b), hi = Math.max(a, b);
        let hit = 0;
        for (const r of S.st.residues) {
          if (ch && r.chainID !== ch) continue;
          if (r.resSeq < lo || r.resSeq > hi) continue;
          out.add(r.rid); hit++;
        }
        if (!hit) bad.push(t);
        continue;
      }
      bad.push(t);
    }
    return { rids: out, bad: bad };
  }

  /** Collapses excluded residues into chain ranges for the chip list. */
  function exclRanges() {
    if (!S.st) return [];
    const list = [];
    for (const r of S.st.residues) if (S.excluded[r.rid]) list.push(r);
    list.sort((a, b) => a.chainID.localeCompare(b.chainID) || a.resSeq - b.resSeq);
    const out = [];
    let run = null;
    for (const r of list) {
      if (run && run.chain === r.chainID && r.resSeq === run.last + 1 && !r.iCode) {
        run.last = r.resSeq; run.rids.push(r.rid);
        continue;
      }
      run = { chain: r.chainID, first: r.resSeq, last: r.resSeq, rids: [r.rid], iCode: r.iCode };
      out.push(run);
    }
    return out.map((r) => ({
      label: r.chain + r.first + (r.iCode || '') + (r.last > r.first ? '&ndash;' + r.last : ''),
      rids: r.rids,
    }));
  }

  function updateExclChips() {
    const box = $('exclChips');
    const ranges = exclRanges();
    if (!ranges.length) {
      box.innerHTML = '<span class="info">none</span>';
      return;
    }
    box.innerHTML = ranges.map((r, i) =>
      '<span class="chip" data-i="' + i + '">' + r.label + '<b title="Put these residues back">&times;</b></span>'
    ).join('');
    box._ranges = ranges;
  }

  function buildChainChips() {
    const box = $('chainChips');
    if (!S.st || S.st.chains.length < 2) { box.innerHTML = ''; return; }
    box.innerHTML = S.st.chains.map((c) =>
      '<button class="chip on" data-chain="' + c.id + '" title="' + c.residues
      + ' residues">' + c.id + '</button>').join('');
  }

  // ------------------------------------------------------------------------------------ events

  function bind() {
    $('structFile').addEventListener('change', (e) => {
      if (e.target.files[0]) handleFile(e.target.files[0], 'struct');
    });
    $('ppmFile').addEventListener('change', (e) => {
      if (e.target.files[0]) handleFile(e.target.files[0], 'ppm');
    });

    // Drag and drop anywhere over the viewport.
    const stage = $('stage');
    let dragDepth = 0;
    window.addEventListener('dragenter', (e) => {
      if (!e.dataTransfer || !Array.from(e.dataTransfer.types).includes('Files')) return;
      e.preventDefault(); dragDepth++; $('dropOverlay').hidden = false;
    });
    window.addEventListener('dragover', (e) => {
      if (!e.dataTransfer) return;
      e.preventDefault(); e.dataTransfer.dropEffect = 'copy';
    });
    window.addEventListener('dragleave', () => {
      dragDepth = Math.max(0, dragDepth - 1);
      if (!dragDepth) $('dropOverlay').hidden = true;
    });
    window.addEventListener('drop', async (e) => {
      e.preventDefault(); dragDepth = 0; $('dropOverlay').hidden = true;
      const files = Array.from((e.dataTransfer && e.dataTransfer.files) || []);
      // A structure first, so a PPM file dropped with it can attach to it.
      const texts = [];
      for (const f of files) texts.push({ f: f, t: await readFile(f) });
      texts.sort((a, b) => (DUM_RE.test(a.t) ? 1 : 0) - (DUM_RE.test(b.t) ? 1 : 0));
      for (const x of texts) await handleFile(x.f);
    });

    const recalc = ['probe', 'points', 'rsaCut', 'patchR', 'incHet', 'incH'];
    for (const id of recalc) $(id).addEventListener('change', () => schedule());
    $('probe').addEventListener('input', () => schedule(400));
    $('rsaCut').addEventListener('input', () => schedule(400));
    $('patchR').addEventListener('input', () => schedule(400));

    $('scaleSel').addEventListener('change', () => {
      $('scaleNote').innerHTML = (Chem.SCALES[$('scaleSel').value] || {}).note || '';
      schedule(0);
    });
    $('aaOnly').addEventListener('change', () => schedule(0));

    for (const id of ['memMode', 'memThick', 'memCentre']) {
      $(id).addEventListener('change', () => { updateMembrane(); schedule(0); });
    }
    $('memHead').addEventListener('change', () => schedule(0));
    $('memOcclude').addEventListener('change', () => schedule(0));
    $('memSkip').addEventListener('change', () => {
      $('memKeepHead').disabled = !$('memSkip').checked || !S.membrane;
      schedule(0);
    });
    $('memKeepHead').addEventListener('change', () => schedule(0));
    $('memSide').addEventListener('change', () => {
      S.ctx.plusLabel = zoneLabels().plusLabel;
      S.ctx.minusLabel = zoneLabels().minusLabel;
      updateSideFilter();
      refillTable(); restyle();
    });
    $('memShow').addEventListener('change', restyle);
    $('memOpacity').addEventListener('input', restyle);

    for (const id of ['repSel', 'surfSel', 'colourSel', 'showExcl', 'showLabel']) {
      $(id).addEventListener('change', restyle);
    }
    $('opacity').addEventListener('input', restyle);
    $('resetBtn').addEventListener('click', () => Viewer.resetView());
    $('pngBtn').addEventListener('click', () => {
      const uri = Viewer.snapshot();
      if (!uri) return;
      const a = document.createElement('a');
      a.href = uri;
      a.download = (S.fileName || 'structure').replace(/\.[^.]+$/, '') + '_epitope.png';
      document.body.appendChild(a); a.click(); a.remove();
    });

    $('showSel').addEventListener('change', refillTable);
    $('sideSel').addEventListener('change', refillTable);
    $('searchBox').addEventListener('input', refillTable);
    $('chainChips').addEventListener('click', (e) => {
      const b = e.target.closest('button[data-chain]');
      if (!b) return;
      b.classList.toggle('on');
      const on = Array.from($('chainChips').querySelectorAll('button.on')).map((x) => x.dataset.chain);
      S.chainFilter = on.length === S.st.chains.length ? null : new Set(on);
      refillTable();
    });
    $('csvBtn').addEventListener('click', () => {
      download((S.fileName || 'structure').replace(/\.[^.]+$/, '') + '_surface_residues.csv',
        Table.csv(), 'text/csv');
    });

    $('exclAdd').addEventListener('click', applyExclText);
    $('exclText').addEventListener('keydown', (e) => { if (e.key === 'Enter') applyExclText(); });
    $('exclClear').addEventListener('click', () => {
      if (!S.st) return;
      S.excluded.fill(0);
      recompute();
    });
    $('exclHet').addEventListener('click', () => {
      if (!S.st) return;
      const rids = [];
      for (const r of S.st.residues) if (!r.isAA) rids.push(r.rid);
      if (!rids.length) { message($('exclMsg'), 'There are no non-amino-acid groups to exclude.'); return; }
      message($('exclMsg'), '');
      setExcluded(rids, true);
    });
    $('exclChips').addEventListener('click', (e) => {
      const b = e.target.closest('b');
      if (!b) return;
      const chip = b.parentElement;
      const ranges = $('exclChips')._ranges || [];
      const r = ranges[parseInt(chip.dataset.i, 10)];
      if (r) setExcluded(r.rids, false);
    });

    $('pickClear').addEventListener('click', () => {
      S.picked.clear(); refillTable(); updateSummary(); restyle();
    });
    $('pickPatch').addEventListener('click', () => {
      if (S.focus === null || !S.byRid[S.focus]) {
        message($('exclMsg'), 'Select a residue in the table first; the patch is taken around it.');
        return;
      }
      const o = analysisOpts();
      const c = S.byRid[S.focus];
      const r2 = o.patchR * o.patchR;
      let n = 0;
      for (const row of S.rows) {
        if (!SASA.isCandidate(row, candOpts())) continue;
        const dx = row.surfPt[0] - c.surfPt[0], dy = row.surfPt[1] - c.surfPt[1],
              dz = row.surfPt[2] - c.surfPt[2];
        if (dx * dx + dy * dy + dz * dz <= r2) { S.picked.add(row.rid); n++; }
      }
      message($('exclMsg'), n + ' residues picked within ' + o.patchR + '&nbsp;&Aring; of '
        + c.label + ' ' + c.chainID + '.');
      refillTable(); updateSummary(); restyle();
    });
    $('copyHot').addEventListener('click', async () => {
      const s = hotspotString();
      if (!s) return;
      const done = () => {
        $('copyHot').textContent = 'Copied';
        setTimeout(() => { $('copyHot').textContent = 'Copy'; }, 1200);
      };
      try {
        await navigator.clipboard.writeText(s);
        done();
      } catch (err) {
        // navigator.clipboard needs a secure context, which a page opened from the file system is
        // not. Fall back to selecting the text and asking the document to copy it.
        $('hotspots').select();
        $('hotspots').setSelectionRange(0, s.length);
        let ok = false;
        try { ok = document.execCommand('copy'); } catch (e2) { ok = false; }
        if (ok) done(); else $('copyHot').textContent = 'Ctrl+C';
      }
    });

    $('savePdb').addEventListener('click', savePDB);
    $('saveJson').addEventListener('click', saveSession);

    // Table height: drag the grip at the top of the dock.
    const grip = $('dockGrip');
    grip.addEventListener('pointerdown', (e) => {
      grip.setPointerCapture(e.pointerId);
      const startY = e.clientY;
      const startH = $('dock').getBoundingClientRect().height;
      const move = (ev) => {
        const h = Math.min(window.innerHeight - 160, Math.max(90, startH - (ev.clientY - startY)));
        document.documentElement.style.setProperty('--dock-h', h + 'px');
        Viewer.resize();
      };
      const up = () => {
        grip.removeEventListener('pointermove', move);
        grip.removeEventListener('pointerup', up);
        Viewer.resize();
      };
      grip.addEventListener('pointermove', move);
      grip.addEventListener('pointerup', up);
    });

    window.addEventListener('resize', () => Viewer.resize());
  }

  function applyExclText() {
    if (!S.st) return;
    const { rids, bad } = parseSelection($('exclText').value);
    if (bad.length) {
      message($('exclMsg'), 'Not found: ' + bad.join(', ')
        + '. Use a chain letter with the residue number, such as A137, or a range such as A20-28.',
        'err');
    } else {
      message($('exclMsg'), '');
    }
    if (rids.size) {
      $('exclText').value = '';
      setExcluded(rids, true);
    }
  }

  function savePDB() {
    if (!S.st) return;
    const mode = $('colourSel').value;
    const scale = Colour.scaleFor(mode, S.ctx);
    const field = scale.field || 'rsa';
    const keep = new Uint8Array(S.st.atoms.length);
    const bv = new Float64Array(S.st.atoms.length);
    for (let i = 0; i < S.st.atoms.length; i++) {
      const a = S.st.atoms[i];
      if (S.excluded[a.rid]) continue;
      keep[i] = 1;
      const row = S.byRid[a.rid];
      const v = row ? row[field] : NaN;
      bv[i] = Number.isFinite(v) ? v : 0;
    }
    const head = '' +
      'REMARK   1 Epitope mapper: B-factor column holds ' + field + '\n' +
      'REMARK   1 probe ' + analysisOpts().probe.toFixed(2) + ' A, '
        + analysisOpts().nPoints + ' points per atom\n' +
      (S.membrane ? 'REMARK   1 bilayer z ' + S.membrane.zLo.toFixed(2) + ' to '
        + S.membrane.zHi.toFixed(2) + '\n' : '') +
      (hotspotString() ? 'REMARK   1 picked ' + hotspotString() + '\n' : '');
    download((S.fileName || 'structure').replace(/\.[^.]+$/, '') + '_' + field + '.pdb',
      head + PDB.toPDBText(S.st, { keep: keep, bValues: bv }), 'chemical/x-pdb');
  }

  function saveSession() {
    if (!S.st) return;
    const o = analysisOpts();
    const excluded = [];
    for (const r of S.st.residues) if (S.excluded[r.rid]) excluded.push(r.tag);
    const session = {
      tool: 'epitope_mapper',
      structure: S.fileName,
      ppm: S.ppmName || null,
      settings: {
        probe: o.probe, points: o.nPoints, surfaceRSA: o.surfaceRSA, patchRadius: o.patchR,
        scale: o.scale, includeHet: o.includeHet, includeHydrogens: o.includeH,
        candidatesAminoAcidsOnly: $('aaOnly').checked, headgroup: o.headgroup,
      },
      electrostatics: {
        model: 'Coulomb, distance-dependent dielectric eps = 4r, Debye screening',
        ionicStrength_M: Elec.DEFAULT_IONIC,
        debyeLength_A: +Elec.debyeLength(Elec.DEFAULT_IONIC).toFixed(2),
        cutoff_A: Elec.DEFAULT_CUTOFF,
        chargeSites: S.elec ? S.elec.sites.n : 0,
        netCharge: S.elec ? S.elec.sites.total : 0,
      },
      membrane: S.membrane ? {
        zLo: S.membrane.zLo, zHi: S.membrane.zHi, thickness: S.membrane.thickness,
        source: S.membrane.source, occludes: $('memOcclude').checked,
        skipBilayerResidues: $('memSkip').checked, plusZ: $('memSide').value,
      } : null,
      excluded: excluded,
      picked: hotspotString().split(',').filter(Boolean),
      pickedSummary: (() => {
        const rows = [];
        for (const rid of S.picked) if (S.byRid[rid]) rows.push(S.byRid[rid]);
        const s = SASA.summarise(rows, o.scale);
        return { residues: s.n, sasa: Math.round(s.sasa), meanHydrophobicity: s.hyd, netCharge: s.charge };
      })(),
    };
    download((S.fileName || 'structure').replace(/\.[^.]+$/, '') + '_epitope_session.json',
      JSON.stringify(session, null, 2), 'application/json');
  }

  // -------------------------------------------------------------------------------------- boot

  function start() {
    const sel = $('scaleSel');
    for (const key in Chem.SCALES) {
      const opt = document.createElement('option');
      opt.value = key;
      opt.textContent = Chem.SCALES[key].label;
      sel.appendChild(opt);
    }
    sel.value = 'kd';
    $('scaleNote').innerHTML = Chem.SCALES.kd.note;

    Table.init({ head: $('tHead'), body: $('tBody'), count: $('tableCount') }, {
      focus: (rid, o) => focusResidue(rid, { centre: o && o.centre, scroll: false }),
      togglePick: togglePick,
      toggleDrop: toggleDrop,
      hover: (rid) => { S.hovered = rid; updateReadout(rid === null ? S.focus : rid); },
    });

    try {
      Viewer.init($('view'), { background: '#0d0d0d' });
      S.viewerReady = true;
      Viewer.on('pick', (rid) => { if (rid !== null) focusResidue(rid, { scroll: true }); });
      Viewer.on('hover', (rid) => { S.hovered = rid; updateReadout(rid === null ? S.focus : rid); });
      Viewer.on('unhover', () => { S.hovered = null; updateReadout(S.focus); });
    } catch (err) {
      message($('fileMsg'), String(err.message || err)
        + ' The table and the exports still work.', 'err');
    }

    bind();
    updateExclChips();
    updateMembrane();
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start);
  else start();

  global.App = S;
})(typeof window !== 'undefined' ? window : globalThis);
