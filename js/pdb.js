/* Structure reader for PDB and mmCIF coordinate files.
 *
 * Returns a flat atom list in file order plus a residue index built on top of it. The atom order is
 * preserved exactly, because the viewer is handed a PDB rendered from the same list and the two are
 * matched by index (see render.js).
 *
 * Membrane proteins: a PPM 3.0 / OPM output file carries DUM dummy atoms that mark the two
 * hydrophobic boundary planes of the bilayer, and its coordinates are rotated so that the membrane
 * normal is the z axis. detectMembrane() recovers the slab from those atoms.
 */
(function (global) {
  'use strict';

  const Chem = global.Chem;

  function pad(s, n) { s = String(s); return s.length >= n ? s : ' '.repeat(n - s.length) + s; }
  function padEnd(s, n) { s = String(s); return s.length >= n ? s : s + ' '.repeat(n - s.length); }
  function fixed(v, w, d) {
    let s = (Number.isFinite(v) ? v : 0).toFixed(d);
    if (s.length > w) s = s.slice(0, w);
    return pad(s, w);
  }

  /** Element from the element column when present, otherwise inferred from the atom name. */
  function elementOf(elemCol, atomName, resName) {
    const e = (elemCol || '').trim().toUpperCase();
    if (e && /^[A-Z][A-Z]?$/.test(e)) return e;
    let n = (atomName || '').trim().toUpperCase().replace(/^[0-9]+/, '');
    if (!n) return 'C';
    // Standard residues use the one-letter element in the first position of the name, so CA is
    // carbon alpha and CD is carbon delta rather than calcium or cadmium. SE is the exception: no
    // amino-acid atom carries that name, and selenomethionine and selenocysteine both use it.
    if (Chem.isAminoAcid(resName)) return n === 'SE' ? 'SE' : n[0];
    const two = n.slice(0, 2);
    if (Chem.TWO_LETTER[two] === 1) return two;
    return n[0];
  }

  /** PDB atom-name column rule: one-letter elements are indented by one unless the name is 4 long. */
  function formatAtomName(name, element) {
    const n = (name || '').trim();
    if (n.length >= 4) return n.slice(0, 4);
    if (element && element.length === 2) return padEnd(n, 4);
    return ' ' + padEnd(n, 3);
  }

  function atomToLine(a, serial) {
    const record = a.isHet ? 'HETATM' : 'ATOM  ';
    return record
      + pad(serial, 5) + ' '
      + formatAtomName(a.name, a.element)
      + padEnd(a.altLoc || ' ', 1)
      + pad(a.resName, 3) + ' '
      + padEnd(a.chainID || 'A', 1)
      + pad(a.resSeq, 4)
      + padEnd(a.iCode || ' ', 1) + '   '
      + fixed(a.x, 8, 3) + fixed(a.y, 8, 3) + fixed(a.z, 8, 3)
      + fixed(a.occ, 6, 2) + fixed(a.b, 6, 2) + '          '
      + pad(a.element || '', 2);
  }

  // ------------------------------------------------------------------------------------- PDB

  function parsePDB(text, name) {
    const lines = text.split(/\r\n|\r|\n/);
    const atoms = [];
    const warnings = [];
    const header = [];           // records worth passing on to the viewer and to exports
    let model = 0, modelsSeen = 0, stopped = false;

    for (let i = 0; i < lines.length; i++) {
      const line = lines[i];
      const rec = line.slice(0, 6);
      if (rec === 'MODEL ') {
        modelsSeen++;
        model = parseInt(line.slice(10, 14), 10) || modelsSeen;
        if (modelsSeen > 1) stopped = true;
        continue;
      }
      if (rec === 'ENDMDL') continue;
      if (rec === 'HELIX ' || rec === 'SHEET ' || rec === 'SSBOND' || rec === 'CRYST1') {
        header.push(line);
        continue;
      }
      if (rec !== 'ATOM  ' && rec !== 'HETATM') continue;
      if (stopped) continue;     // only the first model is analysed

      const x = parseFloat(line.slice(30, 38));
      const y = parseFloat(line.slice(38, 46));
      const z = parseFloat(line.slice(46, 54));
      if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
        warnings.push('Line ' + (i + 1) + ' has unreadable coordinates and was skipped.');
        continue;
      }
      const atomName = line.slice(12, 16).trim();
      const resName = line.slice(17, 20).trim().toUpperCase();
      const occRaw = parseFloat(line.slice(54, 60));
      const bRaw = parseFloat(line.slice(60, 66));
      atoms.push({
        serial: parseInt(line.slice(6, 11), 10) || atoms.length + 1,
        name: atomName,
        altLoc: line.slice(16, 17).trim(),
        resName: resName,
        chainID: line.slice(21, 22).trim() || 'A',
        resSeq: parseInt(line.slice(22, 26), 10),
        iCode: line.slice(26, 27).trim(),
        x: x, y: y, z: z,
        occ: Number.isFinite(occRaw) ? occRaw : 1,
        b: Number.isFinite(bRaw) ? bRaw : 0,
        element: elementOf(line.slice(76, 78), atomName, resName),
        isHet: rec === 'HETATM',
        model: model,
      });
    }
    if (modelsSeen > 1) {
      warnings.push('The file holds ' + modelsSeen + ' models; only the first was read.');
    }
    return { atoms: atoms, warnings: warnings, header: header, format: 'pdb' };
  }

  // ----------------------------------------------------------------------------------- mmCIF

  /** Splits an mmCIF loop row, honouring single, double and semicolon-free quoting. */
  function cifTokens(line) {
    const out = [];
    let i = 0;
    while (i < line.length) {
      const c = line[i];
      if (c === ' ' || c === '\t') { i++; continue; }
      if (c === '"' || c === "'") {
        const end = line.indexOf(c, i + 1);
        if (end < 0) { out.push(line.slice(i + 1)); break; }
        out.push(line.slice(i + 1, end));
        i = end + 1;
      } else {
        let end = i;
        while (end < line.length && line[end] !== ' ' && line[end] !== '\t') end++;
        out.push(line.slice(i, end));
        i = end;
      }
    }
    return out;
  }

  function parseCIF(text, name) {
    const lines = text.split(/\r\n|\r|\n/);
    const atoms = [];
    const warnings = [];
    let i = 0, found = false;

    while (i < lines.length) {
      if (lines[i].trim().toLowerCase() !== 'loop_') { i++; continue; }
      // Collect the tag list of this loop.
      const tags = [];
      let j = i + 1;
      while (j < lines.length && lines[j].trim().startsWith('_')) {
        tags.push(lines[j].trim().split(/\s+/)[0].toLowerCase());
        j++;
      }
      if (!tags.length || !tags[0].startsWith('_atom_site.')) { i = j; continue; }
      found = true;

      const col = {};
      tags.forEach((t, k) => { col[t.replace('_atom_site.', '')] = k; });
      const get = (row, a, b) => {
        let k = col[a];
        if (k === undefined) k = col[b];
        return k === undefined ? '' : (row[k] || '');
      };

      let firstModel = null;
      for (; j < lines.length; j++) {
        const raw = lines[j];
        const t = raw.trim();
        if (!t || t === '#' || t === 'loop_' || t.startsWith('_') || t.startsWith('data_')) break;
        const row = cifTokens(raw);
        if (row.length < tags.length) continue;

        const group = get(row, 'group_pdb');
        const modelNum = parseInt(get(row, 'pdbx_pdb_model_num'), 10) || 1;
        if (firstModel === null) firstModel = modelNum;
        if (modelNum !== firstModel) { warnings.push('Only model ' + firstModel + ' was read.'); break; }

        const alt = get(row, 'label_alt_id');
        const ins = get(row, 'pdbx_pdb_ins_code');
        const resName = (get(row, 'auth_comp_id', 'label_comp_id') || '').toUpperCase();
        const atomName = get(row, 'auth_atom_id', 'label_atom_id').replace(/"/g, '');
        const x = parseFloat(get(row, 'cartn_x')), y = parseFloat(get(row, 'cartn_y')),
              z = parseFloat(get(row, 'cartn_z'));
        if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) continue;
        const occ = parseFloat(get(row, 'occupancy')), b = parseFloat(get(row, 'b_iso_or_equiv'));
        atoms.push({
          serial: parseInt(get(row, 'id'), 10) || atoms.length + 1,
          name: atomName,
          altLoc: (alt === '.' || alt === '?') ? '' : alt,
          resName: resName,
          chainID: get(row, 'auth_asym_id', 'label_asym_id') || 'A',
          resSeq: parseInt(get(row, 'auth_seq_id', 'label_seq_id'), 10),
          iCode: (ins === '.' || ins === '?') ? '' : ins,
          x: x, y: y, z: z,
          occ: Number.isFinite(occ) ? occ : 1,
          b: Number.isFinite(b) ? b : 0,
          element: elementOf(get(row, 'type_symbol'), atomName, resName),
          isHet: group.toUpperCase() === 'HETATM',
          model: modelNum,
        });
      }
      i = j;
    }
    if (!found) throw new Error('"' + name + '" contains no _atom_site records (not a coordinate mmCIF file).');
    return { atoms: atoms, warnings: warnings, header: [], format: 'cif' };
  }

  // ------------------------------------------------------------------- clean-up and indexing

  /** Chain IDs longer than one character cannot be written to a PDB column; map them to letters. */
  function normaliseChains(atoms, warnings) {
    const long = {};
    for (const a of atoms) if (a.chainID.length > 1) long[a.chainID] = true;
    const names = Object.keys(long);
    if (!names.length) return;
    const taken = {};
    for (const a of atoms) if (a.chainID.length === 1) taken[a.chainID] = true;
    const pool = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'.split('');
    const map = {};
    for (const n of names) {
      let pick = n[0];
      if (taken[pick]) pick = pool.find((c) => !taken[c]) || '?';
      taken[pick] = true;
      map[n] = pick;
    }
    for (const a of atoms) if (map[a.chainID]) a.chainID = map[a.chainID];
    warnings.push('Multi-character chain names were shortened: '
      + names.map((n) => n + ' &rarr; ' + map[n]).join(', ') + '.');
  }

  /** Keeps one alternate location per atom: the blank/"A" one, or the highest occupancy. */
  function resolveAltLocs(atoms, warnings) {
    let any = false;
    for (const a of atoms) if (a.altLoc) { any = true; break; }
    if (!any) return atoms;

    const best = new Map();
    for (const a of atoms) {
      const key = a.chainID + '|' + a.resSeq + '|' + a.iCode + '|' + a.resName + '|' + a.name;
      const prev = best.get(key);
      if (!prev) { best.set(key, a); continue; }
      const better = a.occ > prev.occ + 1e-6
        || (Math.abs(a.occ - prev.occ) <= 1e-6 && (a.altLoc === '' || a.altLoc < prev.altLoc));
      if (better) best.set(key, a);
    }
    const keep = new Set(best.values());
    const out = atoms.filter((a) => keep.has(a));
    warnings.push('Alternate locations found: kept ' + out.length + ' of ' + atoms.length
      + ' atoms (highest occupancy per atom).');
    return out;
  }

  /** Builds the residue index. Waters and membrane dummies are routed out of the atom list. */
  function index(parsed, name, opts) {
    const warnings = parsed.warnings.slice();
    normaliseChains(parsed.atoms, warnings);
    let atoms = resolveAltLocs(parsed.atoms, warnings);

    const dummies = [];
    let waters = 0;
    atoms = atoms.filter((a) => {
      if (Chem.isMembraneDummy(a.resName)) { dummies.push(a); return false; }
      if (Chem.isWater(a.resName)) { waters++; return false; }
      return true;
    });
    // A file of nothing but dummy atoms still defines a bilayer, which is a valid thing to load on
    // its own, so the caller can allow it.
    if (!atoms.length && !(opts && opts.allowEmpty && dummies.length)) {
      throw new Error('"' + name + '" has no protein or ligand atoms. '
        + (dummies.length ? 'It holds only DUM membrane dummy atoms; load the structure as well.'
                          : 'Expected a PDB or mmCIF coordinate file.'));
    }

    const residues = [];
    const byKey = new Map();
    let nH = 0, nHet = 0;

    atoms.forEach((a, i) => {
      a.index = i;
      if (a.element === 'H' || a.element === 'D') { a.isH = true; nH++; } else a.isH = false;
      a.radius = Chem.radius(a.element);
      const key = a.chainID + '|' + a.resSeq + '|' + a.iCode;
      let res = byKey.get(key);
      if (!res) {
        const canon = Chem.canonical(a.resName);
        const isAA = Chem.isAminoAcid(a.resName);
        res = {
          rid: residues.length,
          key: key,
          chainID: a.chainID,
          resSeq: a.resSeq,
          iCode: a.iCode,
          resName: a.resName,
          canon: canon,
          one: Chem.one(a.resName),
          isAA: isAA,
          isHet: true,              // flipped to false below if any ATOM record contributes
          atoms: [],
          label: Chem.one(a.resName) + a.resSeq + a.iCode,
          tag: a.chainID + a.resSeq + a.iCode,
          cls: Chem.residueClass(a.resName),
        };
        residues.push(res);
        byKey.set(key, res);
      }
      if (!a.isHet) res.isHet = false;
      a.rid = res.rid;
      res.atoms.push(i);
    });

    for (const r of residues) {
      if (r.isHet) nHet++;
      if (!r.isAA) r.label = r.resName + r.resSeq + r.iCode;
    }

    const chains = [];
    const chainMap = new Map();
    for (const r of residues) {
      let c = chainMap.get(r.chainID);
      if (!c) {
        c = { id: r.chainID, residues: 0, aminoAcids: 0, atoms: 0 };
        chains.push(c); chainMap.set(r.chainID, c);
      }
      c.residues++;
      if (r.isAA) c.aminoAcids++;
      c.atoms += r.atoms.length;
    }

    return {
      name: name,
      format: parsed.format,
      atoms: atoms,
      residues: residues,
      chains: chains,
      dummies: dummies,
      header: parsed.header,
      counts: { atoms: atoms.length, hydrogens: nH, waters: waters, hetResidues: nHet,
                residues: residues.length, dummies: dummies.length },
      warnings: warnings,
    };
  }

  /** Auto-detects the format and returns the indexed structure. opts: { allowEmpty }. */
  function parse(text, name, opts) {
    if (!text || !text.trim()) throw new Error('"' + name + '" is empty.');
    const looksCIF = /^\s*data_/m.test(text) || /_atom_site\./.test(text);
    const looksPDB = /^(ATOM  |HETATM)/m.test(text);
    let parsed;
    if (looksPDB) parsed = parsePDB(text, name);
    else if (looksCIF) parsed = parseCIF(text, name);
    else {
      throw new Error('"' + name + '" is not a coordinate file: no ATOM/HETATM records and no '
        + 'mmCIF atom_site loop were found.');
    }
    if (!parsed.atoms.length) {
      throw new Error('"' + name + '" contains no atoms that could be read.');
    }
    return index(parsed, name, opts);
  }

  /** PDB text for a structure, in atom-list order. opts.bValues overrides the B-factor column. */
  function toPDBText(st, opts) {
    opts = opts || {};
    const keep = opts.keep || null;        // Uint8Array over st.atoms, 1 = write this atom
    const bv = opts.bValues || null;       // Float array over st.atoms
    const out = [];
    if (opts.header !== false) for (const h of st.header) out.push(h);
    let serial = 1, lastChain = null, lastRes = null;
    for (let i = 0; i < st.atoms.length; i++) {
      if (keep && !keep[i]) continue;
      const a = st.atoms[i];
      if (lastChain !== null && a.chainID !== lastChain) {
        out.push('TER   ' + pad(serial++, 5) + '      ' + pad(lastRes.resName, 3) + ' '
          + padEnd(lastChain, 1) + pad(lastRes.resSeq, 4) + padEnd(lastRes.iCode || ' ', 1));
      }
      const rec = bv ? Object.assign({}, a, { b: bv[i] }) : a;
      out.push(atomToLine(rec, serial++));
      lastChain = a.chainID;
      lastRes = st.residues[a.rid];
    }
    if (lastChain !== null) {
      out.push('TER   ' + pad(serial++, 5) + '      ' + pad(lastRes.resName, 3) + ' '
        + padEnd(lastChain, 1) + pad(lastRes.resSeq, 4) + padEnd(lastRes.iCode || ' ', 1));
    }
    out.push('END');
    return out.join('\n') + '\n';
  }

  /** Recovers the bilayer slab from PPM/OPM DUM atoms, which lie on the two boundary planes. */
  function detectMembrane(dummies) {
    if (!dummies || dummies.length < 4) return null;
    const zs = dummies.map((a) => a.z).sort((p, q) => p - q);
    const lo = zs[0], hi = zs[zs.length - 1];
    if (hi - lo < 4) return null;          // a single plane: no slab to define

    // Two-means split on z, seeded at the extremes. PPM writes two flat grids of dummy atoms.
    let cLo = lo, cHi = hi;
    for (let it = 0; it < 25; it++) {
      let sLo = 0, nLo = 0, sHi = 0, nHi = 0;
      for (const z of zs) {
        if (Math.abs(z - cLo) <= Math.abs(z - cHi)) { sLo += z; nLo++; } else { sHi += z; nHi++; }
      }
      const newLo = nLo ? sLo / nLo : cLo, newHi = nHi ? sHi / nHi : cHi;
      if (Math.abs(newLo - cLo) < 1e-6 && Math.abs(newHi - cHi) < 1e-6) break;
      cLo = newLo; cHi = newHi;
    }
    let spreadLo = 0, spreadHi = 0, nLo = 0, nHi = 0;
    for (const z of zs) {
      if (Math.abs(z - cLo) <= Math.abs(z - cHi)) { spreadLo = Math.max(spreadLo, Math.abs(z - cLo)); nLo++; }
      else { spreadHi = Math.max(spreadHi, Math.abs(z - cHi)); nHi++; }
    }
    const xs = dummies.map((a) => a.x), ys = dummies.map((a) => a.y);
    return {
      zLo: cLo,
      zHi: cHi,
      thickness: cHi - cLo,
      centre: (cLo + cHi) / 2,
      normal: 'z',
      planes: [nLo, nHi],
      flat: Math.max(spreadLo, spreadHi) < 1.5,
      extent: {
        xMin: Math.min.apply(null, xs), xMax: Math.max.apply(null, xs),
        yMin: Math.min.apply(null, ys), yMax: Math.max.apply(null, ys),
      },
      source: 'dummy atoms',
    };
  }

  global.PDB = { parse: parse, parsePDB: parsePDB, parseCIF: parseCIF, toPDBText: toPDBText,
                 detectMembrane: detectMembrane, atomToLine: atomToLine, elementOf: elementOf };
})(typeof window !== 'undefined' ? window : globalThis);
