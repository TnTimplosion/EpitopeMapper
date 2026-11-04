# Epitope mapper

Interactive tool for picking **binder epitopes** on a target structure. Load a PDB (or mmCIF), and it
computes the solvent accessible surface area of every residue, lists the exposed ones with their
hydrophobicity and charge, and colours the surface so you can see what a binder would actually meet.
Click a residue in the table and it lights up on the structure; click the structure and the table
follows. For a membrane protein, feed it the **PPM 3.0** output and the bilayer is taken into
account, so lipid-facing residues stop being mistaken for binding sites.

Runs entirely in the browser. Nothing is uploaded, nothing needs installing.

## Start

1. Open `index.html` (or the single-file build `dist/epitope_mapper.html`) in Chrome, Edge or Firefox.
2. Drop your files on the window, or use the buttons in the panel:
   - **Structure**: the target, `.pdb`, `.ent` or `.cif`.
   - **PPM 3.0 output**: optional, for membrane proteins. This is the file the
     [PPM server](https://opm.phar.umich.edu/ppm_server3) returns with `DUM` dummy atoms on the two
     bilayer boundary planes. It also holds the structure rotated into the membrane frame, so loading
     it alone is enough: the app uses those coordinates for everything.
3. Work down the table. Tick **✓** to collect residues into an epitope, **✗** to exclude a residue
   from the analysis entirely.

Mouse: drag = rotate, shift-drag or right-drag = pan, wheel = zoom, click a residue = select it.
In the table: click a row = select, double-click = centre the view on it, click a header = sort,
drag the bar above the table = resize it.

## What you see

| Element | Meaning |
|---|---|
| 3D view | Cartoon plus the selected surface, coloured by the property chosen in **Display › Colour by**. The legend at the bottom left always describes the current colouring |
| Readout (top left) | The residue under the pointer, or the selected one: area, relative area, hydrophobicity, charge, patch area, membrane region |
| Table | One row per residue. Highlighted cells use the same colours as the 3D view. Sorted by relative SASA to start with, so the most exposed residues are at the top |
| **Show** | Which residues to list. *Epitope candidates* is the working set: exposed, not excluded, amino acid, and outside the bilayer |
| **Side** | Only present once a bilayer is defined: lists one face of the membrane at a time. Takes the names you give the faces, so it reads *Extracellular only* rather than *+z side only* once you have said which way round the protein sits |
| Chain chips | Only present for multi-chain structures; click to show or hide a chain's residues |
| Picked epitope | Running totals for the ticked residues, and the hotspot list in the `A30,A33` form that binder pipelines expect |

## The analysis

**Accessible surface.** Shrake-Rupley: test points are spread over a sphere of radius *r*+1.4 Å around
each atom and counted as accessible when no other atom covers them. Radii are Bondi van der Waals
values; the probe radius and the number of test points are both editable. Hydrogens are ignored by
default, waters are always dropped, and ligands, cofactors and glycans occlude the surface unless you
untick that box. 100 points per atom puts the total within ≈0.2 % of the 480-point value; use 480
when you want numbers to quote.

**Relative SASA** is the residue's area as a percentage of the reference maximum for that residue
type (Tien et al. 2013, theoretical). A residue counts as *surface* at RSA ≥ 20 % by default. Non-standard
residues have no reference area, so they fall back to an absolute cut of 10 Å².

**Hydrophobicity** comes from one of four scales. Two of them run in the opposite direction, which the
tool handles for you: colours and patch averages always put hydrophobic at the warm end.

| Scale | Units | Direction |
|---|---|---|
| Kyte-Doolittle (1982) | — | positive = hydrophobic |
| Eisenberg consensus (1984) | — | positive = hydrophobic |
| Wimley-White (1996), water → octanol | kcal/mol | **negative = hydrophobic** |
| Hopp-Woods (1981) hydrophilicity | — | **positive = hydrophilic** |

**Charge** is the formal side-chain charge near pH 7.4: Asp and Glu −1, Lys and Arg +1, His +0.1
(the fraction implied by a side-chain pKa near 6). Chain termini are not charged, because in most
structures they are truncation points rather than real termini.

**Electrostatics.** *Colour by › Electrostatic potential* shades the surface with a Coulombic
potential, the approximation viewers use when a Poisson-Boltzmann solver is not at hand:

```
phi(r) = 332.06 * sum_i  q_i * exp(-d_i / lambda) / (eps(d_i) * d_i)      kcal / (mol e)
```

with a distance-dependent dielectric `eps(d) = 4d` and a Debye factor for the salt in the buffer
(`lambda = 3.04 Å / sqrt(I)`, 7.9 Å at 150 mM; terms past 25 Å are dropped). The charges are the
formal charges of the ionisable groups, spread over the atoms that share them — a carboxylate
carries −0.5 on each oxygen, a guanidinium +0.5 on each of NH1 and NH2, histidine the same +0.1 used
elsewhere — plus any free ions and metal centres in the structure when ligands are switched on.

Two details make the numbers usable rather than singular. The potential for an atom is evaluated on
its *solvent-accessible sphere* rather than at its centre: distances shorter than the atom's radius
plus the probe are clamped to that, which is the closest a water, or a binder atom, can approach. And
only atoms in the SASA calculation carry charge, so excluding a residue removes its charge too — the
map is the one the protein would present without it. Unlike every other colouring, this one is
**per atom**, not per residue, so the surface varies across a single side chain; the table's **ESP**
column is the area-weighted mean over each residue's exposed atoms.

**Patches.** A single residue is not an epitope. For every candidate residue the tool sums the
exposed area, the area-weighted hydrophobicity and the net charge of all candidates within the patch
radius (10 Å by default) of its exposed surface. Sorting by **Patch A** finds the largest contiguous
accessible surfaces; **Patch hyd** and **Patch q** say what kind of surface they are.

## Membrane proteins

When a PPM 3.0 output is loaded, the two planes of `DUM` atoms define the hydrophobic slab, and the
coordinates are already rotated so that the bilayer normal is z. Three things follow from that:

- **The lipid occludes.** A test point that falls inside the slab is not solvent accessible, so
  interfacial residues get the area that is really available to a soluble binder, not the area they
  would have in vacuum. The area lost this way is reported in the **Lipid** column.
- **Residues are assigned a region**: bilayer core, headgroup band (5 Å outside each plane by
  default), or aqueous on the +z or −z side. If you tell the tool which side is extracellular, the
  table says so instead of "+z".
- **Bilayer residues are kept out of the candidate set**, since a soluble binder cannot reach them.
  Headgroup residues can be allowed back in with a tick box.
- **The table can be narrowed to one face**, which is usually the whole point for a membrane
  target: an antibody or minibinder reaches the extracellular side only. The **Side** filter in the
  table toolbar takes each face with its own headgroup band; residues in the bilayer core belong to
  neither and drop out of both. *Export CSV* writes the rows as listed, so it gives you that face
  alone. The picked set is yours and is left as it is.

The bilayer is drawn as two translucent planes at the boundaries, with an opacity slider. They are
flat quads that take no part in depth writing, so the protein surface stays whole where the membrane
crosses it; depth *testing* still applies, so a plane is hidden behind the protein and tints it
where it passes in front. Seen from the default side-on view a flat plane is edge-on and reads as a
line — tilt the view and the sheets open up.

No PPM output? Set the bilayer by hand (thickness and centre along z) — useful if your model is
already membrane-aligned — or leave the membrane off and the whole surface is treated as aqueous.

## Excluding residues

Excluded residues are removed before the surface is computed, so their neighbours gain the area the
excluded residue was covering. It is the same surface the protein would present if those residues
were not in the structure, not a filter applied afterwards. Use it for purification tags, flexible
loops you intend to truncate, a partner subunit you will leave out, or a bound ligand you will not
have.

Type selections into the box (`A12`, `A20-28`, `12-20`, a bare chain letter, or `het` for everything
that is not an amino acid), or click **✗** on a row. Excluded residues appear as chips, stay visible
in the 3D view as a faint wireframe, and are listed in the exports.

## Exports

| Button | What you get |
|---|---|
| **Export CSV** | The rows currently listed, every column, plus the picked and excluded flags |
| **Copy** (hotspot list) | `A30,A33,A35` — the chain+number form used by the `hotspots` field of a BindCraft target, and by RFdiffusion hotspot residues |
| **Save PDB** | The kept atoms, with the current colouring metric written into the B-factor column, so you can colour by it in PyMOL or ChimeraX. Excluded residues are left out. The settings and the picked list are written as REMARK lines |
| **Save session** | Every setting, the exclusions, the picked residues and their totals, as JSON |
| **Save PNG** | The 3D view as it stands |

## Things to know

- The structure is treated as rigid. Side chains do not move, so a residue that is marginally buried
  in this conformation may still be reachable in another.
- The patch numbers are geometric sums over a sphere, not an affinity prediction, and not a
  prediction of immunogenicity. They rank sites by how much accessible surface of what character is
  available; the choice is still yours.
- Charges are formal, not computed from the structure. There is no pKa calculation, so a buried
  carboxylate still reads −1.
- The electrostatic potential is Coulombic, not Poisson-Boltzmann. There is no solute/solvent
  dielectric boundary and no ion-exclusion layer, and the only charges are the ionisable groups, so
  amide and hydroxyl dipoles contribute nothing. It tells you whether a patch is acidic, basic or
  neutral and how strongly; it is not a binding energy. For a publication figure, run APBS.
- Phosphorylated and other modified residues are mapped to their parent amino acid, so a
  phosphoserine carries serine's charge, not the phosphate's.
- The bilayer is an implicit flat slab. Headgroups are a band, not atoms, and the lipid is not
  modelled as a surface: a residue in a groove open to lipid gets no solvent area from it.
- The structure and the dummy atoms must be in the same frame. Load the PPM output itself rather
  than mixing its planes with a differently oriented copy of the structure — the panel warns you if
  no residue ends up inside the slab.
- Secondary structure in the **SS** column is 3Dmol's geometric assignment, not DSSP.
- Only the first model of a multi-model file is read. One alternate location is kept per atom: the
  blank one, or the highest occupancy. mmCIF input reads the `atom_site` loop only, preferring the
  `auth_*` numbering.
- Timings on a laptop: 327 atoms ≈ 20 ms, 5,252 atoms ≈ 150 ms at 100 points per atom. The surface
  mesh, not the area calculation, is what costs time on large structures; set **Surface** to *None*
  while you work and switch it back on at the end. Structures over 20,000 atoms start with the
  surface off.

## Files

```
index.html, css/, js/          the app
  js/chem.js                   radii, reference areas, scales, charges
  js/pdb.js                    PDB and mmCIF reader, PPM dummy-atom detection, PDB writer
  js/sasa.js                   Shrake-Rupley, per-residue metrics, patches
  js/elec.js                   charge sites and the Coulombic surface potential
  js/colour.js                 colour ramps shared by the view, the table and the legend
  js/render.js                 3Dmol.js wrapper: styling, surfaces, membrane, picking
  js/table.js                  the residue table
  js/app.js                    wiring
vendor/3Dmol-min.js            3Dmol.js 2.4.2, unmodified
dist/epitope_mapper.html       single-file build (python tools/build_single_file.py)
tools/serve.py                 local server for the tests, with caching turned off
tools/sasa_reference.py        brute-force SASA reference (numpy only)
tests/                         browser tests (see below)
```

## Tests

`tests/core_test.html` runs from the file system. `tests/accuracy_test.html` reads its fixtures with
`fetch`, so serve the folder first:

```
python tools/serve.py
```

- `http://localhost:8765/tests/core_test.html` — 132 checks: the area of an isolated atom and of an
  overlapping pair against closed-form geometry, the neighbour grid against a brute-force
  implementation, parser behaviour (alternate locations, waters, dummy atoms, element inference,
  round trip), bilayer detection, slab occlusion, exclusion, patches, the Coulombic potential
  against its closed form, group and ion charges, colour scales.
- `http://localhost:8765/tests/accuracy_test.html` — per-residue areas for crambin against
  `tools/sasa_reference.py`, with and without a slab (agreement to 5 × 10⁻⁷ Å²), plus convergence in
  the number of test points.

Regenerate the references with:

```
python tools/sasa_reference.py tests/1crn.pdb -o tests/ref_sasa.csv
python tools/sasa_reference.py tests/1crn.pdb --slab -5 5 -o tests/ref_sasa_slab.csv
```

## Licence

MIT — see [LICENSE](LICENSE). Use it, change it, fork it.

`vendor/3Dmol-min.js` is an unmodified build of [3Dmol.js](https://3dmol.csb.pitt.edu/) 2.4.2,
redistributed under its own BSD-3-Clause licence: see [vendor/3Dmol-LICENSE.txt](vendor/3Dmol-LICENSE.txt).
The test structure `tests/1crn.pdb` is [PDB entry 1CRN](https://www.rcsb.org/structure/1CRN),
redistributed from the public RCSB archive.

## References

- Shrake & Rupley (1973) *J Mol Biol* **79**:351 — the area calculation.
- Tien, Meyer, Sydykova, Spielman & Wilke (2013) *PLoS ONE* **8**:e80635 — reference maximum areas.
- Bondi (1964) *J Phys Chem* **68**:441 — van der Waals radii.
- Kyte & Doolittle (1982) *J Mol Biol* **157**:105; Eisenberg et al. (1984) *J Mol Biol* **179**:125;
  Wimley & White (1996) *Nat Struct Biol* **3**:842; Hopp & Woods (1981) *PNAS* **78**:3824 — scales.
- Lomize, Todd & Lomize (2022) *Protein Sci* **31**:209 — PPM 3.0.
- Pettersen et al. (2021) *Protein Sci* **30**:70 — the Coulombic surface colouring this follows,
  including the distance-dependent dielectric. Baker et al. (2001) *PNAS* **98**:10037 — APBS, for
  when the Poisson-Boltzmann answer is the one you need.
- Rego & Koes (2015) *Bioinformatics* **31**:1322 — 3Dmol.js.
