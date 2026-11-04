/* 3D view: a thin layer over 3Dmol.js.
 *
 * The viewer is handed a PDB rendered from the parsed atom list, in the same order, so every 3Dmol
 * atom can be matched to a parsed atom by position. The match is verified after loading and falls
 * back to a chain/residue/atom-name key if the orders ever disagree.
 *
 * Styling is driven entirely by two properties written onto the 3Dmol atoms:
 *   grp  0 kept, 1 excluded, 2 picked, 3 focused   - selects which representation an atom gets
 *   surf 1 when the atom is part of the SASA calculation - selects what the surface is built from
 * Selecting on a plain property is O(1) per atom in 3Dmol, which keeps restyling fast on large
 * structures, and per-atom colours are written to atom.color (with a colorfunc as a belt-and-braces
 * second route) so the surface, cartoon and table always agree.
 */
(function (global) {
  'use strict';

  const Colour = global.Colour;
  const PICK_COLOUR = 0xfab219;
  const FOCUS_COLOUR = 0x4cd0ff;
  const EXCLUDED_COLOUR = 0x44443f;

  let viewer = null;
  let model = null;
  let st = null;
  let refs = [];                 // refs[i] = 3Dmol atom for parsed atom i
  let keyToRid = new Map();      // "chain|resSeq|iCode" -> residue id, for picking
  let membraneShapes = [];       // the bilayer planes, kept so their materials can be adjusted
  let lastMembrane = null;       // the arguments they were drawn from, for a redraw
  let homeView = null;
  let surfacePending = 0;
  let lastSurfaceKey = '';
  let callbacks = { pick: null, hover: null, unhover: null };

  function available() {
    return typeof global.$3Dmol !== 'undefined' && !!global.$3Dmol.createViewer;
  }

  function init(container, opts) {
    opts = opts || {};
    if (!available()) throw new Error('3Dmol.js did not load, so the 3D view cannot be created.');
    // Surfaces are normally built in web workers created from a blob URL, which a page opened
    // straight from the file system is not allowed to do. The synchronous path gives the same
    // surface on the main thread, so the single-file build works by double-clicking it.
    if (global.location && global.location.protocol === 'file:'
        && typeof global.$3Dmol.setSyncSurface === 'function') {
      global.$3Dmol.setSyncSurface(true);
    }
    viewer = global.$3Dmol.createViewer(container, {
      backgroundColor: opts.background || '#0d0d0d',
      antialias: true,
      id: 'epitope-view',
    });
    viewer.setHoverDuration(120);
    return viewer;
  }

  function clear() {
    if (!viewer) return;
    viewer.removeAllModels();
    viewer.removeAllSurfaces();
    viewer.removeAllShapes();
    viewer.removeAllLabels();
    model = null; st = null; refs = []; keyToRid = new Map(); homeView = null; lastSurfaceKey = '';
    membraneShapes = [];
    viewer.render();
  }

  /** Loads a parsed structure. Returns { matched, secondary } diagnostics. */
  function load(structure) {
    st = structure;
    viewer.removeAllModels();
    viewer.removeAllSurfaces();
    viewer.removeAllShapes();
    viewer.removeAllLabels();

    keyToRid = new Map();
    for (const res of st.residues) keyToRid.set(res.key, res.rid);

    const text = global.PDB.toPDBText(st, {});
    model = viewer.addModel(text, 'pdb', { keepH: true });
    const mAtoms = model.selectedAtoms({});

    let matched = 'index';
    refs = new Array(st.atoms.length).fill(null);
    let ok = mAtoms.length === st.atoms.length;
    if (ok) {
      for (let i = 0; i < mAtoms.length; i++) {
        const a = mAtoms[i], b = st.atoms[i];
        if (a.resi !== b.resSeq || String(a.atom || '').trim() !== b.name) { ok = false; break; }
      }
    }
    if (ok) {
      for (let i = 0; i < mAtoms.length; i++) refs[i] = mAtoms[i];
    } else {
      matched = 'key';
      const byKey = new Map();
      for (const a of mAtoms) {
        const ic = (a.icode || '').trim();
        byKey.set(String(a.chain || '').trim() + '|' + a.resi + '|' + ic + '|' + String(a.atom || '').trim(), a);
      }
      let hits = 0;
      st.atoms.forEach((b, i) => {
        const a = byKey.get(b.chainID + '|' + b.resSeq + '|' + b.iCode + '|' + b.name);
        if (a) { refs[i] = a; hits++; }
      });
      if (!hits) matched = 'failed';
    }

    // Secondary structure, as assigned by 3Dmol from the geometry: majority vote per residue.
    const votes = new Map();
    st.atoms.forEach((b, i) => {
      const a = refs[i];
      if (!a || !a.ss) return;
      let v = votes.get(b.rid);
      if (!v) { v = {}; votes.set(b.rid, v); }
      v[a.ss] = (v[a.ss] || 0) + 1;
    });
    let secondary = 0;
    for (const res of st.residues) {
      const v = votes.get(res.rid);
      if (!v) { res.ss = ''; continue; }
      let best = '', n = 0;
      for (const k in v) if (v[k] > n) { n = v[k]; best = k; }
      res.ss = best === 'h' ? 'H' : best === 's' ? 'E' : 'C';
      if (res.ss !== 'C') secondary++;
    }

    viewer.zoomTo();
    homeView = viewer.getView();
    return { matched: matched, secondary: secondary, atoms: mAtoms.length };
  }

  function loaded() { return !!model; }

  /** Writes grp/surf/colour onto the 3Dmol atoms, then applies the representations. */
  function style(state) {
    if (!viewer || !model) return;
    const rows = state.rows || [];
    const mode = state.mode || 'chain';
    const scale = state.scale || null;
    const ctx = state.ctx || {};
    const picked = state.picked || new Set();
    const focus = state.focus;
    const active = state.active || null;

    const byRid = new Array(st.residues.length);
    for (const r of rows) byRid[r.rid] = r;

    const colourOf = new Int32Array(st.residues.length);
    const grpOf = new Uint8Array(st.residues.length);
    for (const res of st.residues) {
      const row = byRid[res.rid];
      let c = Colour.GREY, g = 0;
      if (row) {
        c = Colour.residueColour(row, mode, scale, ctx);
        if (row.excluded) { g = 1; c = EXCLUDED_COLOUR; }
        if (picked.has(res.rid)) { g = 2; }
        if (focus === res.rid) { g = 3; }
      }
      colourOf[res.rid] = c;
      grpOf[res.rid] = g;
    }

    // Most properties belong to a residue, but a potential belongs to a point in space: when the
    // caller supplies per-atom values, each atom is coloured on its own so the surface varies
    // across a side chain instead of stepping residue by residue.
    const atomValues = state.atomValues || null;
    const atomScale = state.atomScale || null;

    for (let i = 0; i < st.atoms.length; i++) {
      const a = refs[i];
      if (!a) continue;
      const b = st.atoms[i];
      const g = grpOf[b.rid];
      a.grp = g;
      a.surf = active ? (active[i] ? 1 : 0) : 1;
      if (atomValues && atomScale && g !== 1) {
        const v = atomValues[i];
        a.color = Number.isFinite(v) ? atomScale.int(v) : Colour.GREY;
      } else {
        a.color = colourOf[b.rid];
      }
    }

    const colorfunc = (a) => (a && typeof a.color === 'number' ? a.color : Colour.GREY);
    const rep = state.rep || 'cartoon';
    const repStyle = (extra) => {
      const s = {};
      if (rep === 'cartoon') {
        s.cartoon = Object.assign({ colorfunc: colorfunc, thickness: 0.3, arrows: true }, extra);
        s.stick = Object.assign({ radius: 0.08, colorfunc: colorfunc }, extra);
      } else if (rep === 'stick') {
        s.stick = Object.assign({ radius: 0.15, colorfunc: colorfunc }, extra);
      } else if (rep === 'sphere') {
        s.sphere = Object.assign({ scale: 0.9, colorfunc: colorfunc }, extra);
      } else if (rep === 'line') {
        s.line = Object.assign({ linewidth: 1.5, colorfunc: colorfunc }, extra);
      } else if (rep === 'trace') {
        s.cartoon = Object.assign({ colorfunc: colorfunc, style: 'trace', thickness: 0.4 }, extra);
      }
      return s;
    };

    viewer.setStyle({}, {});
    viewer.setStyle({ grp: [0, 2, 3] }, repStyle());

    if (state.showExcluded) {
      viewer.addStyle({ grp: 1 }, { line: { linewidth: 1, color: EXCLUDED_COLOUR } });
    }
    // Picked and focused residues get explicit side-chain geometry so they read through a surface.
    viewer.addStyle({ grp: 2 }, {
      stick: { radius: 0.22, color: PICK_COLOUR },
      sphere: { scale: 0.26, color: PICK_COLOUR },
    });
    viewer.addStyle({ grp: 3 }, {
      stick: { radius: 0.3, color: FOCUS_COLOUR },
      sphere: { scale: 0.34, color: FOCUS_COLOUR },
    });

    // ----- surface
    const surfKind = state.surface || 'none';
    const surfKey = [surfKind, state.opacity, mode, state.scale ? state.scale.domain : '',
                     state.surfaceSelKey || '', state.colourKey || ''].join('|');
    if (surfKind === 'none') {
      viewer.removeAllSurfaces();
      lastSurfaceKey = '';
    } else if (surfKey !== lastSurfaceKey) {
      viewer.removeAllSurfaces();
      lastSurfaceKey = surfKey;
      const T = global.$3Dmol.SurfaceType;
      const type = surfKind === 'VDW' ? T.VDW : surfKind === 'SAS' ? T.SAS : T.MS;
      const sel = { surf: 1, grp: [0, 2, 3] };
      surfacePending++;
      const done = () => {
        surfacePending = Math.max(0, surfacePending - 1);
        viewer.render();
        // 3Dmol blends transparent objects in the order they entered the scene, and rebuilds
        // surfaces before shapes on every pass - so from the second pass on, the planes are laid
        // over the surface. The frame that first receives an off-thread surface mesh is the one
        // exception, and it would paint the surface over the planes, so it is followed by another.
        if (lastMembrane) viewer.render();
      };
      let p;
      try {
        p = viewer.addSurface(type, {
          opacity: state.opacity === undefined ? 0.9 : state.opacity,
          colorfunc: colorfunc,
          clickable: true,
          callback: onSurfaceClick,
        }, sel, sel);
      } catch (e) {
        surfacePending = Math.max(0, surfacePending - 1);
        p = null;
      }
      if (p && typeof p.then === 'function') p.then(done, done); else done();
    }

    // ----- membrane planes
    viewer.removeAllShapes();
    membraneShapes = [];
    lastMembrane = null;
    if (state.membrane && state.showMembrane) {
      drawMembrane(state.membrane, state.membraneOpacity);
    }

    // ----- focus label
    viewer.removeAllLabels();
    if (state.label && focus !== null && focus !== undefined && byRid[focus]) {
      const row = byRid[focus];
      viewer.addLabel(state.label, {
        position: { x: row.surfPt[0], y: row.surfPt[1], z: row.surfPt[2] },
        backgroundColor: '#101010', backgroundOpacity: 0.8,
        fontColor: '#f2f1ea', fontSize: 12, borderThickness: 0.5,
        borderColor: '#4cd0ff', inFront: true, alignment: 'bottomCenter',
      });
    }

    viewer.setClickable({}, true, onAtomClick);
    viewer.setHoverable({}, true, onAtomHover, onAtomUnhover);
    viewer.render();
  }

  /**
   * The bilayer: one flat translucent quad on each boundary plane, with a border so the edges read
   * when the plane is seen face-on.
   *
   * Flat quads rather than thin boxes, because a box blends its two parallel faces and comes out
   * twice as dense as asked for. The shapes are registered in membraneShapes so their materials can
   * be opened up after the first render; see seeThroughMembrane().
   */
  function drawMembrane(mem, opacity) {
    lastMembrane = { mem: mem, opacity: opacity };
    const e = mem.extent;
    const pad = 2;
    const x0 = e.xMin - pad, x1 = e.xMax + pad, y0 = e.yMin - pad, y1 = e.yMax + pad;
    const colour = '#e0c24c';
    const alpha = Number.isFinite(opacity) ? opacity : 0.22;
    membraneShapes = [];

    for (const z of [mem.zLo, mem.zHi]) {
      const corners = [
        { x: x0, y: y0, z: z }, { x: x1, y: y0, z: z },
        { x: x1, y: y1, z: z }, { x: x0, y: y1, z: z },
      ];
      const normals = corners.map(() => ({ x: 0, y: 0, z: 1 }));
      membraneShapes.push(seeThrough(viewer.addCustom({
        vertexArr: corners,
        normalArr: normals,
        faceArr: [0, 1, 2, 0, 2, 3],
        color: colour,
        opacity: alpha,
      })));
      for (let i = 0; i < 4; i++) {
        membraneShapes.push(seeThrough(viewer.addLine({
          start: corners[i], end: corners[(i + 1) % 4], color: colour,
        })));
      }
    }
  }

  /**
   * Lets the protein show through a bilayer plane.
   *
   * 3Dmol draws every transparent object into the depth buffer in one pass before blending any of
   * them, so a translucent plane hides whatever transparent geometry sits behind it: the protein
   * surface vanishes where the membrane crosses it. Clearing depthWrite on the plane's material
   * settles both passes at once, because the depth-only material is cloned from it, and depth
   * *testing* stays on, so geometry in front of a plane still covers it.
   *
   * The catch is that a shape rebuilds its mesh and its material from scratch on every render, so
   * the flag has to be set again each time. Wrapping the shape's own globj does that, and wrapping
   * the instance rather than the prototype leaves the rest of the library alone.
   */
  function seeThrough(shape) {
    if (!shape || typeof shape.globj !== 'function') return shape;
    const build = shape.globj.bind(shape);
    shape.globj = function (group) {
      const out = build(group);
      const obj = shape.shapeObj;
      if (obj && obj.children) {
        for (const child of obj.children) {
          if (child.material) child.material.depthWrite = false;
        }
      }
      return out;
    };
    return shape;
  }

  /**
   * Turns the view side-on to the bilayer, which is how a membrane protein wants to be looked at:
   * the membrane normal is z in a PPM frame, so the slab reads as two horizontal lines.
   */
  function orientMembrane() {
    if (!viewer) return;
    viewer.zoomTo();
    viewer.rotate(90, 'x');
    viewer.render();
    homeView = viewer.getView();
  }

  function ridOfAtom(a) {
    if (!a || !st) return null;
    // 3Dmol atoms carry no parsed index, so go back through the residue key.
    const ic = (a.icode || '').trim();
    const key = String(a.chain || '').trim() + '|' + a.resi + '|' + ic;
    let rid = keyToRid.get(key);
    if (rid === undefined) rid = keyToRid.get(key.slice(0, key.lastIndexOf('|') + 1));
    return rid === undefined ? null : rid;
  }

  function onAtomClick(atom) {
    if (callbacks.pick) callbacks.pick(ridOfAtom(atom), atom);
  }

  function onSurfaceClick(_shape, _vwr, event) {
    // Surfaces report the click position rather than an atom; find the nearest exposed atom.
    if (!callbacks.pick || !st) return;
    const p = event && event.intersect ? event.intersect : null;
    if (!p) return;
    let best = null, bestD = Infinity;
    for (let i = 0; i < st.atoms.length; i++) {
      const a = st.atoms[i];
      const d = (a.x - p.x) ** 2 + (a.y - p.y) ** 2 + (a.z - p.z) ** 2;
      if (d < bestD) { bestD = d; best = a; }
    }
    if (best) callbacks.pick(best.rid, best);
  }

  function onAtomHover(atom, _v, _e, _c) {
    if (callbacks.hover) callbacks.hover(ridOfAtom(atom), atom);
  }

  function onAtomUnhover() {
    if (callbacks.unhover) callbacks.unhover();
  }

  function on(name, fn) { callbacks[name] = fn; }

  function centreOn(rid, zoom) {
    if (!viewer || rid === null || rid === undefined) return;
    const res = st.residues[rid];
    if (!res) return;
    const sel = { chain: res.chainID, resi: res.resSeq };
    try {
      if (zoom) viewer.zoomTo(sel, 400); else viewer.center(sel, 400);
    } catch (e) { /* selection may be empty after a reload */ }
  }

  function resetView() {
    if (!viewer) return;
    if (homeView) viewer.setView(homeView); else viewer.zoomTo();
    viewer.render();
  }

  function resize() { if (viewer) viewer.resize(); }

  function snapshot() { return viewer ? viewer.pngURI() : null; }

  function setBackground(c) { if (viewer) { viewer.setBackgroundColor(c); viewer.render(); } }

  function busy() { return surfacePending > 0; }

  global.Viewer = {
    available, init, clear, load, loaded, style, on, centreOn, resetView, resize, snapshot,
    setBackground, busy, orientMembrane,
    PICK_COLOUR, FOCUS_COLOUR,
    get raw() { return viewer; },
  };
})(typeof window !== 'undefined' ? window : globalThis);
