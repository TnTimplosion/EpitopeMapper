#!/usr/bin/env python3
"""Brute-force SASA reference for the browser implementation.

    python tools/sasa_reference.py tests/1crn.pdb -o tests/ref_sasa.csv

Same definition as js/sasa.js (Shrake-Rupley, Bondi radii, 1.4 A probe, golden-spiral test points)
but a deliberately dumb implementation: every atom pair is tested, with no neighbour grid and no
early exit. Agreement to ~1e-9 therefore checks the grid, the bucket arithmetic and the
last-occluder shortcut rather than the formula.

Pass --slab ZLO ZHI to reproduce the lipid-occlusion mode, where a test point inside the
hydrophobic slab counts as buried.
"""
import argparse
import math
import pathlib
import sys

import numpy as np

# Must match VDW in js/chem.js.
VDW = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'P': 1.80, 'S': 1.80, 'CL': 1.75,
    'BR': 1.85, 'I': 1.98, 'SE': 1.90, 'NA': 2.27, 'MG': 1.73, 'K': 2.75, 'CA': 2.31, 'MN': 1.61,
    'FE': 1.56, 'CO': 1.50, 'NI': 1.63, 'CU': 1.40, 'ZN': 1.39, 'CD': 1.58, 'HG': 1.55, 'B': 1.92,
    'SI': 2.10, 'AS': 1.85, 'D': 1.20,
}
VDW_DEFAULT = 1.70
WATERS = {'HOH', 'WAT', 'DOD', 'H2O', 'TIP3', 'SOL'}
DUMMIES = {'DUM', 'DUMM'}
AA = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
      'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
TWO_LETTER = {'CL', 'BR', 'SE', 'ZN', 'FE', 'MG', 'MN', 'NA', 'CU', 'NI', 'CO', 'CA', 'CD', 'HG'}


def element_of(elem_col, atom_name, res_name):
    e = elem_col.strip().upper()
    if e and e.isalpha():
        return e
    n = atom_name.strip().upper().lstrip('0123456789')
    if not n:
        return 'C'
    if res_name.strip().upper() in AA:
        return n[0]
    if n[:2] in TWO_LETTER:
        return n[:2]
    return n[0]


def read_pdb(path):
    """Atoms in file order, waters and dummies dropped, one alternate location per atom."""
    atoms = []
    stopped = False
    for line in pathlib.Path(path).read_text().splitlines():
        rec = line[:6]
        if rec == 'MODEL ':
            if atoms:
                stopped = True
            continue
        if stopped or rec not in ('ATOM  ', 'HETATM'):
            continue
        res_name = line[17:20].strip().upper()
        if res_name in WATERS or res_name in DUMMIES:
            continue
        name = line[12:16].strip()
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        try:
            occ = float(line[54:60])
        except ValueError:
            occ = 1.0
        atoms.append({
            'name': name,
            'altloc': line[16:17].strip(),
            'res_name': res_name,
            'chain': line[21:22].strip() or 'A',
            'res_seq': int(line[22:26]),
            'icode': line[26:27].strip(),
            'xyz': (x, y, z),
            'occ': occ,
            'element': element_of(line[76:78], name, res_name),
            'het': rec == 'HETATM',
        })

    # Keep one alternate location per atom: blank/first altloc, or the highest occupancy.
    best = {}
    for i, a in enumerate(atoms):
        key = (a['chain'], a['res_seq'], a['icode'], a['res_name'], a['name'])
        prev = best.get(key)
        if prev is None:
            best[key] = i
            continue
        b = atoms[prev]
        better = a['occ'] > b['occ'] + 1e-6 or (
            abs(a['occ'] - b['occ']) <= 1e-6 and (a['altloc'] == '' or a['altloc'] < b['altloc']))
        if better:
            best[key] = i
    keep = set(best.values())
    return [a for i, a in enumerate(atoms) if i in keep]


def sphere_points(n):
    k = np.arange(n)
    y = 1.0 - (2.0 * k + 1.0) / n
    r = np.sqrt(np.maximum(0.0, 1.0 - y * y))
    phi = k * math.pi * (3.0 - math.sqrt(5.0))
    return np.stack([np.cos(phi) * r, y, np.sin(phi) * r], axis=1)


def sasa(atoms, probe=1.4, n_points=100, include_h=False, include_het=True, slab=None):
    xyz = np.array([a['xyz'] for a in atoms], dtype=float)
    rad = np.array([VDW.get(a['element'], VDW_DEFAULT) for a in atoms], dtype=float) + probe
    active = np.array([
        (include_h or a['element'] not in ('H', 'D')) and (include_het or not a['het'])
        for a in atoms])
    pts = sphere_points(n_points)

    out_solvent = np.zeros(len(atoms))
    out_free = np.zeros(len(atoms))
    idx = np.where(active)[0]
    for i in idx:
        ri = rad[i]
        d = np.linalg.norm(xyz[idx] - xyz[i], axis=1)
        near = idx[(d < ri + rad[idx]) & (idx != i)]
        shell = xyz[i] + ri * pts                                  # (n_points, 3)
        if len(near):
            dist = np.linalg.norm(shell[:, None, :] - xyz[near][None, :, :], axis=2)
            free = ~(dist < rad[near][None, :]).any(axis=1)
        else:
            free = np.ones(n_points, dtype=bool)
        unit = 4.0 * math.pi * ri * ri / n_points
        out_free[i] = free.sum() * unit
        if slab is None:
            out_solvent[i] = out_free[i]
        else:
            inside = (shell[:, 2] >= slab[0]) & (shell[:, 2] <= slab[1])
            out_solvent[i] = (free & ~inside).sum() * unit
    return out_solvent, out_free, active


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('pdb')
    ap.add_argument('-o', '--out', help='CSV of per-residue areas')
    ap.add_argument('-p', '--probe', type=float, default=1.4)
    ap.add_argument('-n', '--points', type=int, default=100)
    ap.add_argument('--hydrogens', action='store_true')
    ap.add_argument('--no-het', action='store_true')
    ap.add_argument('--slab', nargs=2, type=float, metavar=('ZLO', 'ZHI'))
    args = ap.parse_args()

    atoms = read_pdb(args.pdb)
    if not atoms:
        sys.exit('No atoms read from %s' % args.pdb)
    solvent, free, active = sasa(atoms, args.probe, args.points, args.hydrogens,
                                 not args.no_het, args.slab)

    residues = {}
    order = []
    for i, a in enumerate(atoms):
        key = (a['chain'], a['res_seq'], a['icode'])
        if key not in residues:
            residues[key] = {'res_name': a['res_name'], 'sasa': 0.0, 'free': 0.0}
            order.append(key)
        if active[i]:
            residues[key]['sasa'] += solvent[i]
            residues[key]['free'] += free[i]

    print('%s: %d atoms (%d in the calculation), %d residues'
          % (pathlib.Path(args.pdb).name, len(atoms), int(active.sum()), len(order)))
    print('probe %.2f A, %d points per atom%s'
          % (args.probe, args.points, '' if args.slab is None else ', slab z %.2f..%.2f' % tuple(args.slab)))
    print('total SASA %.2f A^2  (ignoring the slab: %.2f A^2)' % (solvent.sum(), free.sum()))

    if args.out:
        lines = ['chain,resSeq,icode,resName,sasa_A2,sasa_free_A2']
        for key in order:
            r = residues[key]
            lines.append('%s,%d,%s,%s,%.6f,%.6f'
                         % (key[0], key[1], key[2], r['res_name'], r['sasa'], r['free']))
        out = pathlib.Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text('\n'.join(lines) + '\n')
        print('wrote %s (%d residues)' % (out, len(order)))


if __name__ == '__main__':
    main()
