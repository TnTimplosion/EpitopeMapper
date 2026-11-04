/* Amino-acid and element reference data.
 *
 * Everything here is static lookup data: van der Waals radii for the SASA calculation, reference
 * maximum accessible surface areas for relative-SASA, hydrophobicity scales, and formal side-chain
 * charges. Sources are given next to each table and repeated in the README.
 */
(function (global) {
  'use strict';

  /* Bondi van der Waals radii (J. Phys. Chem. 1964, 68, 441), in Angstroem.
   * Used with a 1.4 A water probe by the Shrake-Rupley routine in sasa.js. */
  const VDW = {
    H: 1.20, C: 1.70, N: 1.55, O: 1.52, F: 1.47, P: 1.80, S: 1.80, CL: 1.75, BR: 1.85, I: 1.98,
    SE: 1.90, NA: 2.27, MG: 1.73, K: 2.75, CA: 2.31, MN: 1.61, FE: 1.56, CO: 1.50, NI: 1.63,
    CU: 1.40, ZN: 1.39, CD: 1.58, HG: 1.55, B: 1.92, SI: 2.10, AS: 1.85, D: 1.20,
  };
  const VDW_DEFAULT = 1.70;

  /* Theoretical maximum accessible surface area per residue, in A^2.
   * Tien et al. (2013) PLoS ONE 8(11): e80635, Table 1 ("theoretical" column).
   * Relative SASA (RSA) = 100 * SASA / MAX_ASA[resName]. */
  const MAX_ASA = {
    ALA: 129, ARG: 274, ASN: 195, ASP: 193, CYS: 167, GLN: 225, GLU: 223, GLY: 104, HIS: 224,
    ILE: 197, LEU: 201, LYS: 236, MET: 224, PHE: 240, PRO: 159, SER: 155, THR: 172, TRP: 285,
    TYR: 263, VAL: 174,
  };

  const THREE_TO_ONE = {
    ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G', HIS: 'H',
    ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S', THR: 'T', TRP: 'W',
    TYR: 'Y', VAL: 'V',
  };

  /* Residue names that mean a standard amino acid under a different protonation state, force field
   * or chemical modification. Mapped so that structures from MD, Amber/CHARMM and selenomethionine
   * crystal forms still get radii, scales and reference areas. */
  const ALIASES = {
    MSE: 'MET', CSO: 'CYS', CME: 'CYS', CSD: 'CYS', OCS: 'CYS', CYX: 'CYS', CYM: 'CYS',
    HIE: 'HIS', HID: 'HIS', HIP: 'HIS', HSD: 'HIS', HSE: 'HIS', HSP: 'HIS',
    ASH: 'ASP', GLH: 'GLU', LYN: 'LYS', ARN: 'ARG', TYM: 'TYR',
    SEP: 'SER', TPO: 'THR', PTR: 'TYR', MLY: 'LYS', ALY: 'LYS', M3L: 'LYS',
    KCX: 'LYS', LLP: 'LYS', PCA: 'GLU', SMC: 'CYS', MHO: 'MET', NEP: 'HIS',
    DAL: 'ALA', DAR: 'ARG', DAS: 'ASP', DCY: 'CYS', DGL: 'GLU', DGN: 'GLN', DHI: 'HIS',
    DIL: 'ILE', DLE: 'LEU', DLY: 'LYS', DPN: 'PHE', DPR: 'PRO', DSN: 'SER', DTH: 'THR',
    DTR: 'TRP', DTY: 'TYR', DVA: 'VAL', MVA: 'VAL', ABA: 'ALA', AIB: 'ALA', NLE: 'LEU',
  };

  const WATERS = { HOH: 1, WAT: 1, DOD: 1, H2O: 1, TIP3: 1, SOL: 1 };

  /* Dummy-atom residue names used by PPM / OPM to mark the membrane boundary planes. */
  const MEMBRANE_DUMMY = { DUM: 1, DUMM: 1 };

  /* Hydrophobicity scales.
   *  hydrophobicHigh: true  -> larger value = more hydrophobic
   *  hydrophobicHigh: false -> larger value = more hydrophilic (scale is inverted for colouring) */
  const SCALES = {
    kd: {
      label: 'Kyte-Doolittle',
      units: '',
      note: 'Kyte &amp; Doolittle (1982) J Mol Biol 157:105. Positive = hydrophobic.',
      hydrophobicHigh: true, lo: -4.5, hi: 4.5,
      values: {
        ALA: 1.8, ARG: -4.5, ASN: -3.5, ASP: -3.5, CYS: 2.5, GLN: -3.5, GLU: -3.5, GLY: -0.4,
        HIS: -3.2, ILE: 4.5, LEU: 3.8, LYS: -3.9, MET: 1.9, PHE: 2.8, PRO: -1.6, SER: -0.8,
        THR: -0.7, TRP: -0.9, TYR: -1.3, VAL: 4.2,
      },
    },
    eisenberg: {
      label: 'Eisenberg consensus',
      units: '',
      note: 'Eisenberg et al. (1984) J Mol Biol 179:125, normalised consensus. Positive = hydrophobic.',
      hydrophobicHigh: true, lo: -2.53, hi: 1.38,
      values: {
        ALA: 0.62, ARG: -2.53, ASN: -0.78, ASP: -0.90, CYS: 0.29, GLN: -0.85, GLU: -0.74,
        GLY: 0.48, HIS: -0.40, ILE: 1.38, LEU: 1.06, LYS: -1.50, MET: 0.64, PHE: 1.19,
        PRO: 0.12, SER: -0.18, THR: -0.05, TRP: 0.81, TYR: 0.26, VAL: 1.08,
      },
    },
    ww: {
      label: 'Wimley-White (octanol)',
      units: 'kcal/mol',
      note: 'Wimley &amp; White (1996) Nat Struct Biol 3:842, water&rarr;octanol transfer free energy '
          + 'of Ac-WL-X-LL peptides, charged side chains. Negative = hydrophobic.',
      hydrophobicHigh: false, lo: -2.09, hi: 3.64,
      values: {
        ALA: 0.50, ARG: 1.81, ASN: 0.85, ASP: 3.64, CYS: -0.02, GLN: 0.77, GLU: 3.63, GLY: 1.15,
        HIS: 2.33, ILE: -1.12, LEU: -1.25, LYS: 2.80, MET: -0.67, PHE: -1.71, PRO: 0.14,
        SER: 0.46, THR: 0.25, TRP: -2.09, TYR: -0.71, VAL: -0.46,
      },
    },
    hw: {
      label: 'Hopp-Woods (hydrophilicity)',
      units: '',
      note: 'Hopp &amp; Woods (1981) PNAS 78:3824, the hydrophilicity scale derived for antigenic '
          + 'determinants. Positive = hydrophilic.',
      hydrophobicHigh: false, lo: -3.4, hi: 3.0,
      values: {
        ALA: -0.5, ARG: 3.0, ASN: 0.2, ASP: 3.0, CYS: -1.0, GLN: 0.2, GLU: 3.0, GLY: 0.0,
        HIS: -0.5, ILE: -1.8, LEU: -1.8, LYS: 3.0, MET: -1.3, PHE: -2.5, PRO: 0.0, SER: 0.3,
        THR: -0.4, TRP: -3.4, TYR: -2.3, VAL: -1.5,
      },
    },
  };

  /* Formal side-chain charge near pH 7.4. Histidine is given the fractional occupancy implied by a
   * side-chain pKa of about 6.0. Chain termini are not charged here: in most deposited or predicted
   * structures they are truncation points rather than real termini. */
  const CHARGE = {
    ASP: -1, GLU: -1, LYS: 1, ARG: 1, HIS: 0.1,
    ALA: 0, ASN: 0, CYS: 0, GLN: 0, GLY: 0, ILE: 0, LEU: 0, MET: 0, PHE: 0, PRO: 0, SER: 0,
    THR: 0, TRP: 0, TYR: 0, VAL: 0,
  };

  /* Coarse chemical class, used for the table and for the "residue type" colouring. */
  const CLASS = {
    ALA: 'aliphatic', GLY: 'aliphatic', ILE: 'aliphatic', LEU: 'aliphatic', MET: 'aliphatic',
    PRO: 'aliphatic', VAL: 'aliphatic',
    PHE: 'aromatic', TRP: 'aromatic', TYR: 'aromatic',
    ASN: 'polar', CYS: 'polar', GLN: 'polar', SER: 'polar', THR: 'polar',
    ASP: 'acidic', GLU: 'acidic',
    ARG: 'basic', HIS: 'basic', LYS: 'basic',
  };

  /* Backbone atom names. Side-chain SASA is everything else. */
  const BACKBONE = { N: 1, CA: 1, C: 1, O: 1, OXT: 1, OT1: 1, OT2: 1, H: 1, HA: 1 };

  /* Two-letter element symbols that can appear left-justified in the atom-name columns of
   * non-standard residues, where the element column is often missing. */
  const TWO_LETTER = {
    CL: 1, BR: 1, SE: 1, ZN: 1, FE: 1, MG: 1, MN: 1, NA: 1, CU: 1, NI: 1, CO: 1, CA: 1, CD: 1,
    HG: 1, SI: 1, AS: 1, LI: 1, BE: 1, AL: 1, CR: 1, MO: 1, AG: 1, AU: 1, PT: 1, PD: 1, SR: 1,
    BA: 1, CS: 1, RB: 1, PB: 1, SN: 1,
  };

  /** Canonical residue name: resolves aliases, upper-cases, trims. */
  function canonical(resName) {
    const r = (resName || '').trim().toUpperCase();
    return ALIASES[r] || r;
  }

  /** True for the 20 standard amino acids (after alias resolution). */
  function isAminoAcid(resName) {
    return THREE_TO_ONE[canonical(resName)] !== undefined;
  }

  function one(resName) {
    return THREE_TO_ONE[canonical(resName)] || 'X';
  }

  function radius(element) {
    const e = (element || '').trim().toUpperCase();
    const r = VDW[e];
    return r === undefined ? VDW_DEFAULT : r;
  }

  function maxASA(resName) {
    return MAX_ASA[canonical(resName)] || 0;
  }

  function hydrophobicity(resName, scaleKey) {
    const s = SCALES[scaleKey] || SCALES.kd;
    const v = s.values[canonical(resName)];
    return v === undefined ? NaN : v;
  }

  function charge(resName) {
    const c = CHARGE[canonical(resName)];
    return c === undefined ? 0 : c;
  }

  function residueClass(resName) {
    return CLASS[canonical(resName)] || 'other';
  }

  function isWater(resName) {
    return WATERS[(resName || '').trim().toUpperCase()] === 1;
  }

  function isMembraneDummy(resName) {
    return MEMBRANE_DUMMY[(resName || '').trim().toUpperCase()] === 1;
  }

  function isBackboneAtom(atomName) {
    return BACKBONE[(atomName || '').trim().toUpperCase()] === 1;
  }

  /** Signed hydrophobicity, flipped so that larger always means more hydrophobic. */
  function hydrophobicSign(scaleKey) {
    const s = SCALES[scaleKey] || SCALES.kd;
    return s.hydrophobicHigh ? 1 : -1;
  }

  global.Chem = {
    VDW, VDW_DEFAULT, MAX_ASA, THREE_TO_ONE, ALIASES, SCALES, CHARGE, CLASS, TWO_LETTER,
    canonical, isAminoAcid, one, radius, maxASA, hydrophobicity, charge, residueClass,
    isWater, isMembraneDummy, isBackboneAtom, hydrophobicSign,
  };
})(typeof window !== 'undefined' ? window : globalThis);
