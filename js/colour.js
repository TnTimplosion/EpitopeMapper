/* Colour ramps shared by the 3D view, the table cells and the legend.
 *
 * Every property the structure can be coloured by resolves to a scale object with a domain, a
 * sampler and the stops needed to draw a colour bar, so the three places that show the same number
 * never disagree about what a colour means.
 */
(function (global) {
  'use strict';

  const Chem = global.Chem;

  const RAMPS = {
    /* hydrophilic -> hydrophobic */
    hydro: ['#2f6fb8', '#8ab6e8', '#e9e6dc', '#e8a46a', '#c2451b'],
    /* buried -> exposed */
    exposure: ['#36362f', '#2c6650', '#199e70', '#7fd3ae', '#dcf3e8'],
    /* negative -> positive */
    charge: ['#d93a2b', '#eaa79f', '#e9e6dc', '#a3b8ea', '#2f55b8'],
    /* low -> high */
    spectrum: ['#3b4cc0', '#7f9fe8', '#e9e6dc', '#f0a860', '#b40426'],
    /* buried in lipid -> solvent exposed */
    lipid: ['#efece4', '#f0d08a', '#e0a93f', '#b8791c'],
  };

  const CHAIN_COLOURS = [
    0x4c92e0, 0x5fb87a, 0xe0894c, 0xb07fd4, 0xd4c04c, 0x4cc0c0, 0xd47f9f, 0x8fa8c0,
    0x9fd44c, 0xd44c4c,
  ];

  const ZONE_COLOURS = {
    'core': 0xe0c24c,       // inside the hydrophobic slab
    'head+': 0xd4a06a,      // headgroup region, +z side
    'head-': 0xb08a6a,      // headgroup region, -z side
    'side+': 0x4c92e0,      // aqueous, +z side
    'side-': 0x8f6fd4,      // aqueous, -z side
    '': 0x8fa8a0,
  };

  const CLASS_COLOURS = {
    aliphatic: 0xd4c88f, aromatic: 0xc2a04c, polar: 0x6fc0b8, acidic: 0xd93a2b,
    basic: 0x2f55b8, other: 0x8a8a84,
  };

  const GREY = 0x55554f;

  function hexToRGB(h) {
    const v = parseInt(h.slice(1), 16);
    return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
  }

  /** Samples a list of hex stops at t in [0,1]. */
  function sampleStops(stops, t) {
    if (!Number.isFinite(t)) return [136, 136, 130];
    t = Math.min(1, Math.max(0, t));
    const seg = (stops.length - 1) * t;
    const i = Math.min(stops.length - 2, Math.floor(seg));
    const f = seg - i;
    const a = hexToRGB(stops[i]), b = hexToRGB(stops[i + 1]);
    return [
      Math.round(a[0] + (b[0] - a[0]) * f),
      Math.round(a[1] + (b[1] - a[1]) * f),
      Math.round(a[2] + (b[2] - a[2]) * f),
    ];
  }

  function toInt(rgb) { return (rgb[0] << 16) | (rgb[1] << 8) | rgb[2]; }
  function toCSS(rgb) { return 'rgb(' + rgb[0] + ',' + rgb[1] + ',' + rgb[2] + ')'; }

  /** Readable text colour for a given background. */
  function inkOn(rgb) {
    const lum = (0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]) / 255;
    return lum > 0.55 ? '#14140f' : '#f2f1ea';
  }

  function linearScale(stops, lo, hi, invert) {
    const span = hi - lo;
    return {
      domain: [lo, hi],
      stops: stops,
      invert: !!invert,
      t: function (v) {
        if (!Number.isFinite(v) || span === 0) return NaN;
        const t = (v - lo) / span;
        return invert ? 1 - t : t;
      },
      rgb: function (v) { return sampleStops(stops, this.t(v)); },
      int: function (v) { return Number.isFinite(v) ? toInt(sampleStops(stops, this.t(v))) : GREY; },
      css: function (v) { return Number.isFinite(v) ? toCSS(sampleStops(stops, this.t(v))) : 'transparent'; },
    };
  }

  /**
   * Builds the scale for a colouring mode.
   * ctx: { scale: hydrophobicity scale key, bRange: [lo,hi], rsaMax, patchRange }
   */
  function scaleFor(mode, ctx) {
    ctx = ctx || {};
    const sc = Chem.SCALES[ctx.scale || 'kd'] || Chem.SCALES.kd;
    switch (mode) {
      case 'hydro': {
        const s = linearScale(RAMPS.hydro, sc.lo, sc.hi, !sc.hydrophobicHigh);
        s.label = sc.label + (sc.units ? ' (' + sc.units + ')' : '');
        s.loLabel = sc.hydrophobicHigh ? 'hydrophilic' : 'hydrophobic';
        s.hiLabel = sc.hydrophobicHigh ? 'hydrophobic' : 'hydrophilic';
        s.field = 'hyd';
        return s;
      }
      case 'patchHyd': {
        const s = linearScale(RAMPS.hydro, sc.lo, sc.hi, !sc.hydrophobicHigh);
        s.label = 'Patch mean, ' + sc.label;
        s.loLabel = sc.hydrophobicHigh ? 'hydrophilic' : 'hydrophobic';
        s.hiLabel = sc.hydrophobicHigh ? 'hydrophobic' : 'hydrophilic';
        s.field = 'patchHyd';
        return s;
      }
      case 'rsa': {
        const s = linearScale(RAMPS.exposure, 0, Number.isFinite(ctx.rsaMax) ? ctx.rsaMax : 100, false);
        s.label = 'Relative SASA (%)';
        s.loLabel = 'buried'; s.hiLabel = 'exposed'; s.field = 'rsa';
        return s;
      }
      case 'sasa': {
        const s = linearScale(RAMPS.exposure, 0, Number.isFinite(ctx.sasaMax) ? ctx.sasaMax : 150, false);
        s.label = 'SASA (&Aring;&sup2;)';
        s.loLabel = 'buried'; s.hiLabel = 'exposed'; s.field = 'sasa';
        return s;
      }
      case 'charge': {
        const s = linearScale(RAMPS.charge, -1, 1, false);
        s.label = 'Formal charge at pH 7.4';
        s.loLabel = 'negative'; s.hiLabel = 'positive'; s.field = 'charge';
        return s;
      }
      case 'esp': {
        // Symmetric about zero, so the midpoint of the ramp is always neutral.
        const lim = Number.isFinite(ctx.espMax) ? Math.max(0.5, ctx.espMax) : 5;
        const s = linearScale(RAMPS.charge, -lim, lim, false);
        s.label = 'Coulombic potential (kcal/mol&middot;e)';
        s.loLabel = 'acidic'; s.hiLabel = 'basic';
        s.field = 'esp';
        s.perAtom = true;
        return s;
      }
      case 'patchCharge': {
        const lim = Number.isFinite(ctx.patchChargeMax) ? Math.max(1, ctx.patchChargeMax) : 4;
        const s = linearScale(RAMPS.charge, -lim, lim, false);
        s.label = 'Patch net charge';
        s.loLabel = 'negative'; s.hiLabel = 'positive';
        s.field = 'patchCharge';
        return s;
      }
      case 'patchSASA': {
        const s = linearScale(RAMPS.exposure, 0,
          Number.isFinite(ctx.patchSASAMax) ? ctx.patchSASAMax : 600, false);
        s.label = 'Patch SASA (&Aring;&sup2;)';
        s.loLabel = 'small'; s.hiLabel = 'large'; s.field = 'patchSASA';
        return s;
      }
      case 'lipid': {
        const s = linearScale(RAMPS.lipid, 0,
          Number.isFinite(ctx.lipidMax) ? Math.max(10, ctx.lipidMax) : 120, false);
        s.label = 'Area buried in the bilayer (&Aring;&sup2;)';
        s.loLabel = 'none'; s.hiLabel = 'lipid-facing'; s.field = 'lipid';
        return s;
      }
      case 'bfactor': {
        const r = ctx.bRange || [0, 100];
        const s = linearScale(RAMPS.spectrum, r[0], r[1], false);
        s.label = ctx.bLabel || 'B-factor';
        s.loLabel = ''; s.hiLabel = ''; s.field = 'b';
        return s;
      }
      default: {
        // Categorical modes have no numeric domain; the legend is drawn from keys().
        const s = {
          domain: null, field: null, label: mode === 'zone' ? 'Membrane region'
            : mode === 'cls' ? 'Residue type' : 'Chain',
          int: function () { return GREY; },
          css: function () { return 'transparent'; },
          categorical: true,
        };
        return s;
      }
    }
  }

  /** Legend entries for the categorical modes. */
  function categories(mode, ctx) {
    if (mode === 'zone') {
      const out = [
        { key: 'core', label: 'bilayer core' },
        { key: 'head+', label: 'headgroups, +z' },
        { key: 'head-', label: 'headgroups, &minus;z' },
        { key: 'side+', label: (ctx && ctx.plusLabel) || 'aqueous, +z' },
        { key: 'side-', label: (ctx && ctx.minusLabel) || 'aqueous, &minus;z' },
      ];
      return out.map((e) => Object.assign(e, { colour: ZONE_COLOURS[e.key] }));
    }
    if (mode === 'cls') {
      return ['aliphatic', 'aromatic', 'polar', 'acidic', 'basic'].map((k) => ({
        key: k, label: k, colour: CLASS_COLOURS[k],
      }));
    }
    if (mode === 'chain') {
      const chains = (ctx && ctx.chains) || [];
      return chains.map((c, i) => ({
        key: c.id, label: 'chain ' + c.id, colour: CHAIN_COLOURS[i % CHAIN_COLOURS.length],
      }));
    }
    return [];
  }

  /** Colour of one residue row under a mode, as a 3Dmol integer. */
  function residueColour(row, mode, scale, ctx) {
    switch (mode) {
      case 'chain': {
        const idx = ctx && ctx.chainIndex ? (ctx.chainIndex[row.chainID] || 0) : 0;
        return CHAIN_COLOURS[idx % CHAIN_COLOURS.length];
      }
      case 'cls': return CLASS_COLOURS[row.cls] || CLASS_COLOURS.other;
      case 'zone': return ZONE_COLOURS[row.zone] !== undefined ? ZONE_COLOURS[row.zone] : ZONE_COLOURS[''];
      default: {
        if (!scale || !scale.field) return GREY;
        const v = row[scale.field];
        return scale.int(v);
      }
    }
  }

  global.Colour = {
    RAMPS, CHAIN_COLOURS, ZONE_COLOURS, CLASS_COLOURS, GREY,
    sampleStops, hexToRGB, toInt, toCSS, inkOn, linearScale, scaleFor, categories, residueColour,
    /** CSS gradient string for a colour bar. */
    gradient: function (stops) { return 'linear-gradient(90deg, ' + stops.join(', ') + ')'; },
  };
})(typeof window !== 'undefined' ? window : globalThis);
