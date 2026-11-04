/* Solvent accessible surface area, and the per-residue metrics built from it.
 *
 * SASA uses the numerical Shrake-Rupley method (J Mol Biol 1973, 79:351): test points are spread
 * over a sphere of radius r+probe around each atom and counted as accessible when no other atom's
 * expanded sphere contains them. Neighbours come from a uniform grid, so cost is linear in atoms.
 *
 * Two things are counted in the same pass:
 *   sasa      - accessible to water, with the lipid slab (when defined) treated as occluding
 *   sasaFree  - accessible ignoring the lipid slab
 * Their difference is the area a residue buries in the bilayer, which is what distinguishes a
 * lipid-facing residue from one a soluble binder can actually reach.
 *
 * Excluded residues are dropped from both the probe set and the occluder set, so the result is the
 * surface the protein would present if those residues were not there.
 */
(function (global) {
  'use strict';

  const Chem = global.Chem;
  const PROBE = 1.4;

  /** Evenly spread unit vectors (Fibonacci / golden-spiral sphere). */
  function spherePoints(n) {
    const p = new Float64Array(n * 3);
    const ga = Math.PI * (3 - Math.sqrt(5));
    for (let k = 0; k < n; k++) {
      const y = 1 - (2 * k + 1) / n;
      const r = Math.sqrt(Math.max(0, 1 - y * y));
      const phi = k * ga;
      p[3 * k] = Math.cos(phi) * r;
      p[3 * k + 1] = y;
      p[3 * k + 2] = Math.sin(phi) * r;
    }
    return p;
  }

  /** Uniform bucket grid over the active atoms, for neighbour queries within `cell`. */
  function buildGrid(x, y, z, active, n, cell) {
    let xMin = Infinity, yMin = Infinity, zMin = Infinity;
    let xMax = -Infinity, yMax = -Infinity, zMax = -Infinity, count = 0;
    for (let i = 0; i < n; i++) {
      if (!active[i]) continue;
      count++;
      if (x[i] < xMin) xMin = x[i]; if (x[i] > xMax) xMax = x[i];
      if (y[i] < yMin) yMin = y[i]; if (y[i] > yMax) yMax = y[i];
      if (z[i] < zMin) zMin = z[i]; if (z[i] > zMax) zMax = z[i];
    }
    if (!count) return null;
    const inv = 1 / cell;
    const nx = Math.max(1, Math.floor((xMax - xMin) * inv) + 1);
    const ny = Math.max(1, Math.floor((yMax - yMin) * inv) + 1);
    const nz = Math.max(1, Math.floor((zMax - zMin) * inv) + 1);
    const nCells = nx * ny * nz;
    const counts = new Int32Array(nCells + 1);
    const cellOf = new Int32Array(n).fill(-1);

    for (let i = 0; i < n; i++) {
      if (!active[i]) continue;
      const ix = Math.min(nx - 1, Math.floor((x[i] - xMin) * inv));
      const iy = Math.min(ny - 1, Math.floor((y[i] - yMin) * inv));
      const iz = Math.min(nz - 1, Math.floor((z[i] - zMin) * inv));
      const c = (iz * ny + iy) * nx + ix;
      cellOf[i] = c;
      counts[c + 1]++;
    }
    for (let c = 0; c < nCells; c++) counts[c + 1] += counts[c];
    const start = counts;                      // start[c] .. start[c+1] after the fill below
    const items = new Int32Array(count);
    const cursor = new Int32Array(nCells);
    for (let i = 0; i < n; i++) {
      const c = cellOf[i];
      if (c < 0) continue;
      items[start[c] + cursor[c]++] = i;
    }
    return { xMin, yMin, zMin, nx, ny, nz, inv, start, items };
  }

  /**
   * Starts a chunked SASA calculation. Drive it with step(budgetMs) until it returns true.
   *
   * opts: probe, nPoints, includeH, includeHet, excluded (Uint8Array over residues),
   *       slab ({zLo, zHi} or null)
   */
  function start(st, opts) {
    opts = opts || {};
    const probe = Number.isFinite(opts.probe) ? opts.probe : PROBE;
    const nPoints = Math.max(12, Math.round(opts.nPoints || 100));
    const slab = opts.slab || null;
    const excluded = opts.excluded || null;
    const atoms = st.atoms;
    const n = atoms.length;

    const x = new Float64Array(n), y = new Float64Array(n), z = new Float64Array(n);
    const rr = new Float64Array(n);                 // r + probe
    const active = new Uint8Array(n);
    let maxR = 0, nActive = 0;

    for (let i = 0; i < n; i++) {
      const a = atoms[i];
      x[i] = a.x; y[i] = a.y; z[i] = a.z;
      rr[i] = a.radius + probe;
      let on = true;
      if (a.isH && !opts.includeH) on = false;
      if (a.isHet && !opts.includeHet) on = false;
      if (excluded && excluded[a.rid]) on = false;
      active[i] = on ? 1 : 0;
      if (on) { nActive++; if (rr[i] > maxR) maxR = rr[i]; }
    }

    const sasa = new Float64Array(n);
    const sasaFree = new Float64Array(n);
    const points = spherePoints(nPoints);
    const grid = nActive ? buildGrid(x, y, z, active, n, 2 * maxR) : null;
    const neigh = new Int32Array(Math.min(n, 4096));
    const job = {
      i: 0, n: n, nActive: nActive, nPoints: nPoints, probe: probe,
      sasa: sasa, sasaFree: sasaFree, active: active, done: !grid,
      progress: 0,
    };
    if (!grid) return job;

    const { xMin, yMin, zMin, nx, ny, nz, inv, start: cStart, items } = grid;
    let nbuf = neigh;

    job.step = function (budgetMs) {
      const t0 = (typeof performance !== 'undefined' ? performance.now() : Date.now());
      const deadline = t0 + (budgetMs || 25);
      while (job.i < n) {
        const i = job.i;
        if (!active[i]) { job.i++; continue; }

        const Ri = rr[i], xi = x[i], yi = y[i], zi = z[i];
        // Neighbours: atoms whose expanded sphere can reach this one's test shell.
        let m = 0;
        const ix = Math.min(nx - 1, Math.max(0, Math.floor((xi - xMin) * inv)));
        const iy = Math.min(ny - 1, Math.max(0, Math.floor((yi - yMin) * inv)));
        const iz = Math.min(nz - 1, Math.max(0, Math.floor((zi - zMin) * inv)));
        for (let dz = -1; dz <= 1; dz++) {
          const cz = iz + dz; if (cz < 0 || cz >= nz) continue;
          for (let dy = -1; dy <= 1; dy++) {
            const cy = iy + dy; if (cy < 0 || cy >= ny) continue;
            for (let dx = -1; dx <= 1; dx++) {
              const cx = ix + dx; if (cx < 0 || cx >= nx) continue;
              const c = (cz * ny + cy) * nx + cx;
              for (let k = cStart[c]; k < cStart[c + 1]; k++) {
                const j = items[k];
                if (j === i) continue;
                const dxx = x[j] - xi, dyy = y[j] - yi, dzz = z[j] - zi;
                const cut = Ri + rr[j];
                if (dxx * dxx + dyy * dyy + dzz * dzz >= cut * cut) continue;
                if (m >= nbuf.length) {
                  const bigger = new Int32Array(nbuf.length * 2);
                  bigger.set(nbuf); nbuf = bigger;
                }
                nbuf[m++] = j;
              }
            }
          }
        }

        let accFree = 0, accSolvent = 0, hint = 0;
        for (let p = 0; p < nPoints; p++) {
          const px = xi + Ri * points[3 * p];
          const py = yi + Ri * points[3 * p + 1];
          const pz = zi + Ri * points[3 * p + 2];
          let blocked = false;
          // Start at the neighbour that blocked the previous point: points are spatially coherent.
          for (let q = 0; q < m; q++) {
            const j = nbuf[(q + hint) % m];
            const dxx = x[j] - px, dyy = y[j] - py, dzz = z[j] - pz;
            const rj = rr[j];
            if (dxx * dxx + dyy * dyy + dzz * dzz < rj * rj) {
              blocked = true; hint = (q + hint) % m; break;
            }
          }
          if (blocked) continue;
          accFree++;
          if (slab && pz >= slab.zLo && pz <= slab.zHi) continue;    // inside the bilayer core
          accSolvent++;
        }

        const unit = 4 * Math.PI * Ri * Ri / nPoints;
        sasa[i] = accSolvent * unit;
        sasaFree[i] = accFree * unit;
        job.i++;
        if ((job.i & 255) === 0 && (typeof performance !== 'undefined' ? performance.now() : Date.now()) > deadline) break;
      }
      job.progress = n ? job.i / n : 1;
      job.done = job.i >= n;
      return job.done;
    };
    return job;
  }

  /** Runs a job to completion on the current thread. Used by the tests. */
  function compute(st, opts) {
    const job = start(st, opts);
    if (job.step) while (!job.step(1e9));
    return job;
  }

  // ------------------------------------------------------------------- per-residue aggregation

  /** Which part of the bilayer a residue sits in. Returns '' when no membrane is defined. */
  function membraneZone(zLo, zHi, headgroup, zMin, zMax, zMid) {
    if (!Number.isFinite(zLo)) return '';
    if (zMax >= zLo && zMin <= zHi) return 'core';
    if (zMid > zHi && zMid <= zHi + headgroup) return 'head+';
    if (zMid < zLo && zMid >= zLo - headgroup) return 'head-';
    return zMid > zHi ? 'side+' : 'side-';
  }

  /**
   * Which side of the bilayer a zone sits on: 'plus' or 'minus' for the two aqueous faces and
   * their headgroup bands, 'core' for the interior, '' when there is no membrane. Which face is
   * extracellular is the user's to say, so the names stay geometric here.
   */
  function zoneSide(zone) {
    if (zone === 'side+' || zone === 'head+') return 'plus';
    if (zone === 'side-' || zone === 'head-') return 'minus';
    if (zone === 'core') return 'core';
    return '';
  }

  /**
   * Builds the per-residue table rows.
   *
   * opts: scale, surfaceRSA, excluded, slab, headgroup, includeH, includeHet
   */
  function residueMetrics(st, job, opts) {
    opts = opts || {};
    const scale = opts.scale || 'kd';
    const excluded = opts.excluded || null;
    const slab = opts.slab || null;
    const headgroup = Number.isFinite(opts.headgroup) ? opts.headgroup : 5;
    const surfaceRSA = Number.isFinite(opts.surfaceRSA) ? opts.surfaceRSA : 20;
    const sasa = job.sasa, sasaFree = job.sasaFree, active = job.active;
    const rows = [];

    for (const res of st.residues) {
      let total = 0, free = 0, side = 0, nActive = 0, nExposed = 0;
      let cx = 0, cy = 0, cz = 0, sx = 0, sy = 0, sz = 0, sw = 0, bSum = 0, bN = 0;
      let zMin = Infinity, zMax = -Infinity;

      for (const ai of res.atoms) {
        const a = st.atoms[ai];
        if (a.z < zMin) zMin = a.z;
        if (a.z > zMax) zMax = a.z;
        if (!a.isH) { bSum += a.b; bN++; }
        if (!active[ai]) continue;
        nActive++;
        total += sasa[ai];
        free += sasaFree[ai];
        if (!Chem.isBackboneAtom(a.name)) side += sasa[ai];
        cx += a.x; cy += a.y; cz += a.z;
        if (sasa[ai] > 0) {
          nExposed++;
          sx += a.x * sasa[ai]; sy += a.y * sasa[ai]; sz += a.z * sasa[ai]; sw += sasa[ai];
        }
      }

      const maxA = Chem.maxASA(res.resName);
      const rsa = maxA ? 100 * total / maxA : NaN;
      const zMid = nActive ? cz / nActive : (zMin + zMax) / 2;
      const zone = slab ? membraneZone(slab.zLo, slab.zHi, headgroup, zMin, zMax, zMid) : '';
      const h = Chem.hydrophobicity(res.resName, scale);

      rows.push({
        rid: res.rid,
        res: res,
        chainID: res.chainID,
        resSeq: res.resSeq,
        iCode: res.iCode,
        resName: res.resName,
        canon: res.canon,
        one: res.one,
        label: res.label,
        tag: res.tag,
        cls: res.cls,
        isAA: res.isAA,
        isHet: res.isHet,
        excluded: !!(excluded && excluded[res.rid]),
        nAtoms: res.atoms.length,
        sasa: total,
        sasaFree: free,
        lipid: Math.max(0, free - total),
        scSASA: side,
        rsa: rsa,
        hyd: h,
        charge: Chem.charge(res.resName),
        b: bN ? bSum / bN : NaN,
        zone: zone,
        centroid: nActive ? [cx / nActive, cy / nActive, cz / nActive] : [0, 0, 0],
        surfPt: sw > 0 ? [sx / sw, sy / sw, sz / sw]
                       : (nActive ? [cx / nActive, cy / nActive, cz / nActive] : [0, 0, 0]),
        nExposedAtoms: nExposed,
        // Surface residue: enough relative area to be reachable by a binder. Non-standard residues
        // have no reference area, so they fall back to an absolute cut.
        surface: Number.isFinite(rsa) ? rsa >= surfaceRSA : total >= 10,
        // Patch fields, filled by patchMetrics().
        patchSASA: 0, patchHyd: NaN, patchCharge: 0, patchN: 0,
      });
    }
    return rows;
  }

  /** True when a row can be part of an epitope: exposed, kept, and not shielded by the bilayer. */
  function isCandidate(row, opts) {
    if (row.excluded || !row.surface) return false;
    if (opts && opts.proteinOnly && !row.isAA) return false;
    if (opts && opts.skipMembrane && (row.zone === 'core' || row.zone === 'head+' || row.zone === 'head-')) {
      return opts.keepHeadgroup && row.zone !== 'core';
    }
    return true;
  }

  /**
   * Neighbourhood sums over the candidate set: for each candidate residue, the exposed area,
   * area-weighted hydrophobicity and net charge of every candidate within `radius` of it.
   * This is what makes a patch rather than a single residue comparable between sites.
   */
  function patchMetrics(rows, radius, opts) {
    const cand = rows.filter((r) => isCandidate(r, opts));
    const n = cand.length;
    for (const r of rows) { r.patchSASA = 0; r.patchHyd = NaN; r.patchCharge = 0; r.patchN = 0; }
    if (!n) return;

    const px = new Float64Array(n), py = new Float64Array(n), pz = new Float64Array(n);
    cand.forEach((r, i) => { px[i] = r.surfPt[0]; py[i] = r.surfPt[1]; pz[i] = r.surfPt[2]; });
    const all = new Uint8Array(n).fill(1);
    const grid = buildGrid(px, py, pz, all, n, Math.max(2, radius));
    const r2 = radius * radius;

    for (let i = 0; i < n; i++) {
      let sumA = 0, sumAH = 0, sumQ = 0, count = 0;
      const ix = Math.min(grid.nx - 1, Math.max(0, Math.floor((px[i] - grid.xMin) * grid.inv)));
      const iy = Math.min(grid.ny - 1, Math.max(0, Math.floor((py[i] - grid.yMin) * grid.inv)));
      const iz = Math.min(grid.nz - 1, Math.max(0, Math.floor((pz[i] - grid.zMin) * grid.inv)));
      for (let dz = -1; dz <= 1; dz++) {
        const cz = iz + dz; if (cz < 0 || cz >= grid.nz) continue;
        for (let dy = -1; dy <= 1; dy++) {
          const cy = iy + dy; if (cy < 0 || cy >= grid.ny) continue;
          for (let dx = -1; dx <= 1; dx++) {
            const cx = ix + dx; if (cx < 0 || cx >= grid.nx) continue;
            const c = (cz * grid.ny + cy) * grid.nx + cx;
            for (let k = grid.start[c]; k < grid.start[c + 1]; k++) {
              const j = grid.items[k];
              const ddx = px[j] - px[i], ddy = py[j] - py[i], ddz = pz[j] - pz[i];
              if (ddx * ddx + ddy * ddy + ddz * ddz > r2) continue;
              const rj = cand[j];
              sumA += rj.sasa;
              if (Number.isFinite(rj.hyd)) sumAH += rj.sasa * rj.hyd;
              sumQ += rj.charge;
              count++;
            }
          }
        }
      }
      const r = cand[i];
      r.patchSASA = sumA;
      r.patchHyd = sumA > 0 ? sumAH / sumA : NaN;
      r.patchCharge = sumQ;
      r.patchN = count;
    }
  }

  /** Summary of an arbitrary set of rows: what the picked epitope looks like as a whole. */
  function summarise(rows, scaleKey) {
    let sasa = 0, charge = 0, wh = 0, n = 0, nAA = 0;
    let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity, zMin = Infinity, zMax = -Infinity;
    for (const r of rows) {
      n++;
      if (r.isAA) nAA++;
      sasa += r.sasa;
      charge += r.charge;
      if (Number.isFinite(r.hyd)) wh += r.sasa * r.hyd;
      const p = r.surfPt;
      if (p[0] < xMin) xMin = p[0]; if (p[0] > xMax) xMax = p[0];
      if (p[1] < yMin) yMin = p[1]; if (p[1] > yMax) yMax = p[1];
      if (p[2] < zMin) zMin = p[2]; if (p[2] > zMax) zMax = p[2];
    }
    const span = n ? Math.sqrt((xMax - xMin) ** 2 + (yMax - yMin) ** 2 + (zMax - zMin) ** 2) : 0;
    return {
      n: n, nAA: nAA, sasa: sasa, charge: charge,
      hyd: sasa > 0 ? wh / sasa : NaN,
      span: span, scale: scaleKey,
    };
  }

  global.SASA = {
    PROBE, start, compute, spherePoints, residueMetrics, patchMetrics, isCandidate, summarise,
    membraneZone, zoneSide, buildGrid,
  };
})(typeof window !== 'undefined' ? window : globalThis);
