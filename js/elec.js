/* Electrostatics at the protein surface.
 *
 * A Coulombic potential with a distance-dependent dielectric, which is the approximation molecular
 * viewers use for surface charge when a Poisson-Boltzmann solver is not at hand:
 *
 *     phi(r) = 332.06 * sum_i  q_i * exp(-d_i / lambda) / (eps(d_i) * d_i)      kcal / (mol e)
 *
 * with eps(d) = 4d, which damps long-range terms the way bulk water would, and a Debye factor for
 * the salt in the buffer (lambda = 3.04 A / sqrt(I), 7.9 A at 150 mM). Charges are the formal
 * charges of the ionisable groups, spread over the atoms that share them, so a carboxylate carries
 * -0.5 on each oxygen rather than -1 on the carbon.
 *
 * Only atoms taking part in the SASA calculation carry charge, so excluding a residue removes its
 * charge as well: the potential is the one the protein would present without it.
 *
 * What this is not: it is not a Poisson-Boltzmann calculation. There is no solute/solvent dielectric
 * boundary, no ion exclusion layer, and no pKa shifts - a buried aspartate still carries -1. It
 * ranks the character of a surface patch; it does not give binding energies.
 */
(function (global) {
  'use strict';

  const Chem = global.Chem;
  const COULOMB = 332.0637;          // kcal A / (mol e^2)
  const DEBYE_1M = 3.04;             // Debye length in A at 1 M, 1:1 electrolyte, 25 C
  const DEFAULT_IONIC = 0.15;        // M, physiological buffer
  const DEFAULT_CUTOFF = 25;         // A; exp(-d/7.9)/d^2 is negligible beyond this

  /* Formal charge of each ionisable group, spread over the atoms that carry it. Histidine gets the
   * same fractional occupancy used elsewhere in the app (side-chain pKa near 6, so +0.1 at pH 7.4),
   * split between its two ring nitrogens. */
  const SITES = {
    ASP: { OD1: -0.5, OD2: -0.5 },
    GLU: { OE1: -0.5, OE2: -0.5 },
    LYS: { NZ: 1.0 },
    ARG: { NH1: 0.5, NH2: 0.5 },
    HIS: { ND1: 0.05, NE2: 0.05 },
  };

  /* Free ions and metal centres, keyed by residue name. These sit in the structure as HETATM
   * records and dominate the local potential, so they count when ligands are part of the surface. */
  const IONS = {
    ZN: 2, MG: 2, CA: 2, MN: 2, FE: 2, FE2: 2, FE3: 3, CU: 2, CU1: 1, NI: 2, CO: 2, CD: 2,
    NA: 1, K: 1, LI: 1, RB: 1, CS: 1, SR: 2, BA: 2, AL: 3, HG: 2, PB: 2,
    CL: -1, BR: -1, IOD: -1, F: -1, PO4: -3, SO4: -2, NO3: -1, ACT: -1, FMT: -1, CIT: -3,
  };

  /** Debye screening length in Angstroem for an ionic strength in mol/L. 0 means no screening. */
  function debyeLength(ionicStrength) {
    if (!(ionicStrength > 0)) return Infinity;
    return DEBYE_1M / Math.sqrt(ionicStrength);
  }

  /**
   * Charge sites of a structure. Only atoms flagged active (in the SASA calculation, so not
   * hydrogens, not ligands when those are switched off, and not excluded residues) contribute.
   */
  function sites(st, active) {
    const x = [], y = [], z = [], q = [];
    for (let i = 0; i < st.atoms.length; i++) {
      if (active && !active[i]) continue;
      const a = st.atoms[i];
      const res = Chem.canonical(a.resName);
      let charge = 0;
      const group = SITES[res];
      if (group) {
        const c = group[a.name.toUpperCase()];
        if (c !== undefined) charge = c;
      } else if (IONS[res] !== undefined && st.residues[a.rid].atoms.length === 1) {
        // A one-atom het group named after an ion really is that ion.
        charge = IONS[res];
      }
      if (!charge) continue;
      x.push(a.x); y.push(a.y); z.push(a.z); q.push(charge);
    }
    return {
      n: q.length,
      x: Float64Array.from(x), y: Float64Array.from(y), z: Float64Array.from(z),
      q: Float64Array.from(q),
      total: q.reduce((s, v) => s + v, 0),
    };
  }

  /**
   * Potential at one point, in kcal/(mol e).
   *
   * opts.dmin is the closest approach allowed to a charge. Evaluating at an atom centre, the atom's
   * own charge is singular and a neighbour in the same side chain is unrealistically close, so the
   * caller passes the atom's solvent-accessible radius: the potential is then the one on that
   * atom's accessible sphere, which is where a water, or a binder atom, actually sits.
   */
  function potentialAt(px, py, pz, s, opts) {
    opts = opts || {};
    const lambda = debyeLength(opts.ionicStrength === undefined ? DEFAULT_IONIC : opts.ionicStrength);
    const cutoff = opts.cutoff || DEFAULT_CUTOFF;
    const cut2 = cutoff * cutoff;
    const dmin = opts.dmin > 0 ? opts.dmin : 0;
    let phi = 0;
    for (let i = 0; i < s.n; i++) {
      const dx = s.x[i] - px, dy = s.y[i] - py, dz = s.z[i] - pz;
      const d2 = dx * dx + dy * dy + dz * dz;
      if (d2 > cut2) continue;
      const d = Math.max(Math.sqrt(d2), dmin);
      if (d < 1e-6) continue;
      // eps(d) = 4d, so the 1/(eps*d) term is 1/(4 d^2).
      phi += s.q[i] * Math.exp(-d / lambda) / (4 * d * d);
    }
    return COULOMB * phi;
  }

  /**
   * Potential at every atom that has any surface at all, as a Float64Array over st.atoms.
   * Buried and inactive atoms are left as NaN so they colour as "no value" rather than zero.
   */
  function compute(st, job, opts) {
    opts = opts || {};
    const n = st.atoms.length;
    const phi = new Float64Array(n).fill(NaN);
    const active = job ? job.active : null;
    const s = sites(st, active);
    const evaluated = (i) => !(active && !active[i])
      && !(job && job.sasaFree && !(job.sasaFree[i] > 0));
    if (!s.n) {
      // Nothing ionisable in the structure: the potential is zero, not unknown.
      for (let i = 0; i < n; i++) if (evaluated(i)) phi[i] = 0;
      return { phi: phi, sites: s };
    }

    const lambda = debyeLength(opts.ionicStrength === undefined ? DEFAULT_IONIC : opts.ionicStrength);
    const cutoff = opts.cutoff || DEFAULT_CUTOFF;
    const cut2 = cutoff * cutoff;
    const probe = Number.isFinite(opts.probe) ? opts.probe : 1.4;
    const all = new Uint8Array(s.n).fill(1);
    const grid = global.SASA.buildGrid(s.x, s.y, s.z, all, s.n, cutoff);

    for (let i = 0; i < n; i++) {
      // Only atoms with some exposure can colour a surface; the rest would cost time for nothing.
      if (!evaluated(i)) continue;
      const a = st.atoms[i];
      // Evaluated on this atom's solvent-accessible sphere, not at its centre: see potentialAt().
      const dmin = a.radius + probe;
      let sum = 0;
      const ix = Math.min(grid.nx - 1, Math.max(0, Math.floor((a.x - grid.xMin) * grid.inv)));
      const iy = Math.min(grid.ny - 1, Math.max(0, Math.floor((a.y - grid.yMin) * grid.inv)));
      const iz = Math.min(grid.nz - 1, Math.max(0, Math.floor((a.z - grid.zMin) * grid.inv)));
      for (let dz = -1; dz <= 1; dz++) {
        const cz = iz + dz; if (cz < 0 || cz >= grid.nz) continue;
        for (let dy = -1; dy <= 1; dy++) {
          const cy = iy + dy; if (cy < 0 || cy >= grid.ny) continue;
          for (let dx = -1; dx <= 1; dx++) {
            const cx = ix + dx; if (cx < 0 || cx >= grid.nx) continue;
            const c = (cz * grid.ny + cy) * grid.nx + cx;
            for (let k = grid.start[c]; k < grid.start[c + 1]; k++) {
              const j = grid.items[k];
              const ddx = s.x[j] - a.x, ddy = s.y[j] - a.y, ddz = s.z[j] - a.z;
              const d2 = ddx * ddx + ddy * ddy + ddz * ddz;
              if (d2 > cut2) continue;
              const d = Math.max(Math.sqrt(d2), dmin);
              sum += s.q[j] * Math.exp(-d / lambda) / (4 * d * d);
            }
          }
        }
      }
      phi[i] = COULOMB * sum;
    }
    return { phi: phi, sites: s };
  }

  /**
   * Per-residue potential: the area-weighted mean over the residue's exposed atoms, which is the
   * number the table shows. Residues with no exposed atom get NaN.
   */
  function perResidue(st, job, phi) {
    const out = new Float64Array(st.residues.length).fill(NaN);
    if (!job) return out;
    for (const res of st.residues) {
      let w = 0, sum = 0;
      for (const ai of res.atoms) {
        const area = job.sasa[ai];
        if (!(area > 0) || !Number.isFinite(phi[ai])) continue;
        sum += phi[ai] * area;
        w += area;
      }
      if (w > 0) out[res.rid] = sum / w;
    }
    return out;
  }

  global.Elec = {
    COULOMB, SITES, IONS, DEFAULT_IONIC, DEFAULT_CUTOFF,
    debyeLength, sites, potentialAt, compute, perResidue,
  };
})(typeof window !== 'undefined' ? window : globalThis);
