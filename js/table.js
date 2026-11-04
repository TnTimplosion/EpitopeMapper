/* The surface-residue table.
 *
 * Rows are the per-residue metrics from sasa.js. The table owns sorting, filtering and the two
 * per-row controls that feed back into the analysis: "pick" marks a residue as part of the epitope,
 * "drop" excludes it so the surface is recomputed as if the residue were not in the structure.
 * Clicking anywhere else on a row focuses it in the 3D view.
 */
(function (global) {
  'use strict';

  const Colour = global.Colour;
  const Chem = global.Chem;
  const MAX_ROWS = 2000;

  const COLUMNS = [
    { key: 'pick', label: '&#10003;', title: 'Mark as part of the epitope', cls: 'c-tick', sort: false },
    { key: 'drop', label: '&#10005;', title: 'Exclude this residue from the analysis', cls: 'c-tick', sort: false },
    { key: 'label', label: 'Residue', title: 'One-letter code, number and insertion code', cls: 'c-res' },
    { key: 'chainID', label: 'Ch', title: 'Chain' },
    { key: 'ss', label: 'SS', title: 'Secondary structure assigned from the geometry: H helix, E strand, C coil' },
    { key: 'cls', label: 'Type', title: 'Side-chain chemistry' },
    { key: 'sasa', label: 'SASA', title: 'Solvent accessible surface area, A^2' },
    { key: 'rsa', label: 'RSA %', title: 'SASA as a percentage of the reference maximum for this residue type' },
    { key: 'scSASA', label: 'Side ch.', title: 'Side-chain part of the SASA, A^2' },
    { key: 'hyd', label: 'Hyd', title: 'Hydrophobicity on the selected scale' },
    { key: 'charge', label: 'q', title: 'Formal side-chain charge at pH 7.4' },
    { key: 'esp', label: 'ESP', title: 'Coulombic potential over the exposed atoms of this residue, kcal/mol.e' },
    { key: 'patchSASA', label: 'Patch A', title: 'Summed SASA of the surface residues within the patch radius' },
    { key: 'patchHyd', label: 'Patch hyd', title: 'Area-weighted mean hydrophobicity of the patch' },
    { key: 'patchCharge', label: 'Patch q', title: 'Net charge of the patch' },
    { key: 'lipid', label: 'Lipid', title: 'Area buried in the bilayer, A^2' },
    { key: 'zone', label: 'Region', title: 'Position relative to the membrane' },
  ];

  let els = null;
  let handlers = {};
  let rows = [];
  let view = [];
  let state = {};
  let sortKey = 'rsa';
  let sortDir = -1;
  let focusRid = null;

  function fmt(v, d) {
    if (v === null || v === undefined || (typeof v === 'number' && !Number.isFinite(v))) return '&ndash;';
    return v.toFixed(d);
  }

  function fmtCharge(q) {
    if (!Number.isFinite(q)) return '&ndash;';
    if (q === 0) return '0';
    const s = q > 0 ? '+' : '&minus;';
    const a = Math.abs(q);
    return s + (Number.isInteger(a) ? a.toFixed(0) : a.toFixed(1));
  }

  function zoneText(zone) {
    if (!zone) return '&ndash;';
    const labels = state.zoneLabels || {};
    return labels[zone] || zone;
  }

  function init(elements, h) {
    els = elements;
    handlers = h || {};

    els.head.addEventListener('click', (e) => {
      const th = e.target.closest('th');
      if (!th || !th.dataset.key) return;
      const col = COLUMNS.find((c) => c.key === th.dataset.key);
      if (!col || col.sort === false) return;
      if (sortKey === col.key) sortDir = -sortDir;
      else { sortKey = col.key; sortDir = (col.key === 'label' || col.key === 'chainID' || col.key === 'cls' || col.key === 'zone' || col.key === 'ss') ? 1 : -1; }
      render();
    });

    els.body.addEventListener('click', (e) => {
      const tr = e.target.closest('tr');
      if (!tr || !tr.dataset.rid) return;
      const rid = parseInt(tr.dataset.rid, 10);
      const cell = e.target.closest('td');
      const kind = cell && cell.dataset.act;
      if (kind === 'pick') {
        if (handlers.togglePick) handlers.togglePick(rid, e.shiftKey);
        return;
      }
      if (kind === 'drop') {
        if (handlers.toggleDrop) handlers.toggleDrop(rid, e.shiftKey);
        return;
      }
      if (handlers.focus) handlers.focus(rid, { centre: e.detail > 1 });
    });

    els.body.addEventListener('mouseover', (e) => {
      const tr = e.target.closest('tr');
      if (tr && tr.dataset.rid && handlers.hover) handlers.hover(parseInt(tr.dataset.rid, 10));
    });
    els.body.addEventListener('mouseleave', () => { if (handlers.hover) handlers.hover(null); });
  }

  function matchesSearch(r, q) {
    if (!q) return true;
    const needle = q.toLowerCase();
    return (r.tag.toLowerCase().indexOf(needle) >= 0)
      || (r.label.toLowerCase().indexOf(needle) >= 0)
      || (r.resName.toLowerCase().indexOf(needle) >= 0)
      || (String(r.resSeq).indexOf(needle) >= 0);
  }

  function filter() {
    const show = state.show || 'candidates';
    const q = (state.search || '').trim();
    const chains = state.chains || null;      // Set of chain ids, null = all
    const side = state.side || 'all';         // 'plus' / 'minus' face of the bilayer, or all
    const picked = state.picked || new Set();
    const candOpts = state.candOpts || {};

    return rows.filter((r) => {
      if (chains && !chains.has(r.chainID)) return false;
      // Residues in the bilayer core belong to neither face, so a side filter leaves them out.
      if (side !== 'all' && global.SASA.zoneSide(r.zone) !== side) return false;
      if (!matchesSearch(r, q)) return false;
      switch (show) {
        case 'candidates': return global.SASA.isCandidate(r, candOpts);
        case 'surface': return r.surface && !r.excluded;
        case 'picked': return picked.has(r.rid);
        case 'excluded': return r.excluded;
        case 'buried': return !r.surface && !r.excluded;
        default: return true;
      }
    });
  }

  function compare(a, b) {
    let va = a[sortKey], vb = b[sortKey];
    if (sortKey === 'label') {
      const c = a.chainID.localeCompare(b.chainID);
      if (c) return c * sortDir;
      if (a.resSeq !== b.resSeq) return (a.resSeq - b.resSeq) * sortDir;
      return a.iCode.localeCompare(b.iCode) * sortDir;
    }
    if (typeof va === 'string' || typeof vb === 'string') {
      return String(va || '').localeCompare(String(vb || '')) * sortDir;
    }
    const na = Number.isFinite(va), nb = Number.isFinite(vb);
    if (!na && !nb) return 0;
    if (!na) return 1;            // missing values always sort last
    if (!nb) return -1;
    if (va === vb) {
      const c = a.chainID.localeCompare(b.chainID);
      return c || (a.resSeq - b.resSeq);
    }
    return (va - vb) * sortDir;
  }

  function header() {
    const parts = [];
    for (const c of COLUMNS) {
      const active = sortKey === c.key;
      const arrow = active ? (sortDir > 0 ? ' &uarr;' : ' &darr;') : '';
      parts.push('<th data-key="' + c.key + '" title="' + c.title + '"'
        + (c.cls ? ' class="' + c.cls + (active ? ' sorted' : '') + '"' : (active ? ' class="sorted"' : ''))
        + (c.sort === false ? ' data-nosort="1"' : '')
        + '>' + c.label + arrow + '</th>');
    }
    return '<tr>' + parts.join('') + '</tr>';
  }

  function rowHTML(r) {
    const picked = (state.picked || new Set()).has(r.rid);
    const hydScale = state.hydScale, rsaScale = state.rsaScale, qScale = state.qScale;
    const cls = [];
    if (r.excluded) cls.push('dropped');
    if (picked) cls.push('picked');
    if (r.rid === focusRid) cls.push('focused');
    if (!r.surface) cls.push('buried');

    const tint = (scale, v) => (scale && Number.isFinite(v))
      ? ' style="background:' + scale.css(v) + ';color:' + Colour.inkOn(scale.rgb(v)) + '"' : '';

    return '<tr data-rid="' + r.rid + '" class="' + cls.join(' ') + '">'
      + '<td class="c-tick" data-act="pick"><i class="box' + (picked ? ' on' : '') + '"></i></td>'
      + '<td class="c-tick" data-act="drop"><i class="box' + (r.excluded ? ' off' : '') + '"></i></td>'
      + '<td class="c-res">' + r.label + (r.isHet ? '<sup title="not a standard amino acid">het</sup>' : '') + '</td>'
      + '<td>' + r.chainID + '</td>'
      + '<td class="c-ss ss-' + (r.res.ss || 'C') + '">' + (r.res.ss || '&ndash;') + '</td>'
      + '<td class="c-cls">' + r.cls.slice(0, 4) + '</td>'
      + '<td>' + fmt(r.sasa, 1) + '</td>'
      + '<td' + tint(rsaScale, r.rsa) + '>' + fmt(r.rsa, 1) + '</td>'
      + '<td>' + fmt(r.scSASA, 1) + '</td>'
      + '<td' + tint(hydScale, r.hyd) + '>' + fmt(r.hyd, 2) + '</td>'
      + '<td' + tint(qScale, r.charge) + '>' + fmtCharge(r.charge) + '</td>'
      + '<td' + tint(state.espScale, r.esp) + '>' + fmt(r.esp, 1) + '</td>'
      + '<td>' + fmt(r.patchSASA, 0) + '</td>'
      + '<td' + tint(hydScale, r.patchHyd) + '>' + fmt(r.patchHyd, 2) + '</td>'
      + '<td>' + fmtCharge(r.patchCharge) + '</td>'
      + '<td>' + (r.lipid > 0.5 ? fmt(r.lipid, 0) : '&ndash;') + '</td>'
      + '<td class="c-zone">' + zoneText(r.zone) + '</td>'
      + '</tr>';
  }

  function render() {
    if (!els) return;
    view = filter().sort(compare);
    const shown = view.slice(0, MAX_ROWS);
    els.head.innerHTML = header();
    els.body.innerHTML = shown.map(rowHTML).join('');
    if (els.count) {
      const total = rows.length;
      const picked = (state.picked || new Set()).size;
      const dropped = rows.reduce((n, r) => n + (r.excluded ? 1 : 0), 0);
      els.count.innerHTML = view.length + ' of ' + total + ' residues'
        + (view.length > MAX_ROWS ? ', showing the first ' + MAX_ROWS : '')
        + (picked ? ' &middot; <b>' + picked + ' picked</b>' : '')
        + (dropped ? ' &middot; ' + dropped + ' excluded' : '');
    }
  }

  function setRows(newRows, newState) {
    rows = newRows || [];
    state = newState || {};
    render();
  }

  function update(newState) {
    Object.assign(state, newState || {});
    render();
  }

  function setFocus(rid, scroll) {
    focusRid = rid;
    if (!els) return;
    for (const tr of els.body.querySelectorAll('tr.focused')) tr.classList.remove('focused');
    if (rid === null || rid === undefined) return;
    const tr = els.body.querySelector('tr[data-rid="' + rid + '"]');
    if (tr) {
      tr.classList.add('focused');
      if (scroll) tr.scrollIntoView({ block: 'nearest' });
    } else if (scroll) {
      // The focused residue is filtered out of the current view: say so rather than silently nothing.
      if (handlers.offView) handlers.offView(rid);
    }
  }

  function visibleRows() { return view; }

  function csv() {
    const head = ['chain', 'resSeq', 'insCode', 'resName', 'oneLetter', 'secondary', 'class',
      'sasa_A2', 'rsa_percent', 'sidechain_sasa_A2', 'hydrophobicity', 'hydrophobicity_scale',
      'charge_pH7.4', 'coulombic_potential_kcal_per_mol_e',
      'patch_sasa_A2', 'patch_hydrophobicity', 'patch_charge', 'patch_residues',
      'lipid_buried_A2', 'membrane_region', 'surface', 'picked', 'excluded'];
    const picked = state.picked || new Set();
    const scaleKey = state.scale || 'kd';
    const lines = [head.join(',')];
    for (const r of view) {
      lines.push([
        r.chainID, r.resSeq, r.iCode || '', r.resName, r.one, r.res.ss || '', r.cls,
        r.sasa.toFixed(2), Number.isFinite(r.rsa) ? r.rsa.toFixed(2) : '', r.scSASA.toFixed(2),
        Number.isFinite(r.hyd) ? r.hyd.toFixed(3) : '', scaleKey,
        r.charge, Number.isFinite(r.esp) ? r.esp.toFixed(3) : '', r.patchSASA.toFixed(2),
        Number.isFinite(r.patchHyd) ? r.patchHyd.toFixed(3) : '', r.patchCharge.toFixed(1), r.patchN,
        r.lipid.toFixed(2), r.zone || '', r.surface ? 1 : 0, picked.has(r.rid) ? 1 : 0,
        r.excluded ? 1 : 0,
      ].join(','));
    }
    return lines.join('\n') + '\n';
  }

  global.Table = { init, setRows, update, render, setFocus, csv, visibleRows, COLUMNS, MAX_ROWS };
})(typeof window !== 'undefined' ? window : globalThis);
