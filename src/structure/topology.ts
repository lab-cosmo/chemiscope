/* eslint-disable @typescript-eslint/no-dynamic-delete */
/* eslint-disable @typescript-eslint/no-array-delete */
/* eslint-disable @typescript-eslint/no-unsafe-member-access */
/* eslint-disable @typescript-eslint/no-unsafe-assignment */
/* eslint-disable @typescript-eslint/no-non-null-assertion */
/// Bond assignment from position, taken directly from 3Dmol.js & adapted to
/// TypeScript

/// 3Dmol.js  is licensed under a BSD-3-Clause license.
/// Copyright (c) 2014, University of Pittsburgh and contributors

import { AtomSpec } from '3dmol';

// Covalent radii
const COVALENT_RADII: Record<string, number> = {
    H: 0.37,
    He: 0.32,
    Li: 1.34,
    Be: 0.9,
    B: 0.82,
    C: 0.77,
    N: 0.75,
    O: 0.73,
    F: 0.71,
    Ne: 0.69,
    Na: 1.54,
    Mg: 1.3,
    Al: 1.18,
    Si: 1.11,
    P: 1.06,
    S: 1.02,
    Cl: 0.99,
    Ar: 0.97,
    K: 1.96,
    Ca: 1.74,
    Sc: 1.44,
    Ti: 1.56,
    V: 1.25,
    Mn: 1.39,
    Fe: 1.25,
    Co: 1.26,
    Ni: 1.21,
    Cu: 1.38,
    Zn: 1.31,
    Ga: 1.26,
    Ge: 1.22,
    Se: 1.16,
    Br: 1.14,
    Kr: 1.1,
    Rb: 2.11,
    Sr: 1.92,
    Y: 1.62,
    Zr: 1.48,
    Nb: 1.37,
    Mo: 1.45,
    Tc: 1.56,
    Ru: 1.26,
    Rh: 1.35,
    Pd: 1.31,
    Ag: 1.53,
    Cd: 1.48,
    In: 1.44,
    Sn: 1.41,
    Sb: 1.38,
    Te: 1.35,
    I: 1.33,
    Xe: 1.3,
    Cs: 2.25,
    Ba: 1.98,
    Lu: 1.6,
    Hf: 1.5,
    Ta: 1.38,
    W: 1.46,
    Re: 1.59,
    Os: 1.44,
    Ir: 1.37,
    Pt: 1.28,
    Au: 1.44,
    Hg: 1.49,
    Tl: 1.48,
    Pb: 1.47,
    Bi: 1.46,
    Rn: 1.45,

    // None of the bottom row or any of the Lanthanides have bond lengths
};

function bondLength(elem: string): number {
    return COVALENT_RADII[elem] || 1.6;
}

// return true if atom1 and atom2 are probably bonded to each other
// based on distance alone
function areConnected(atom1: AtomSpec, atom2: AtomSpec) {
    if (atom1.elem === 'X' || atom2.elem === 'X') return false;
    let maxsq = bondLength(atom1.elem!) + bondLength(atom2.elem!);
    maxsq += 0.25; // fudge factor, especially important for md frames, also see 1i3d
    maxsq *= maxsq;

    let xdiff = atom1.x! - atom2.x!;
    xdiff *= xdiff;
    if (xdiff > maxsq) return false;
    let ydiff = atom1.y! - atom2.y!;
    ydiff *= ydiff;
    if (ydiff > maxsq) return false;
    let zdiff = atom1.z! - atom2.z!;
    zdiff *= zdiff;
    if (zdiff > maxsq) return false;

    const distSquared = xdiff + ydiff + zdiff;

    if (isNaN(distSquared)) {
        return false;
    } else if (distSquared < 0.5) {
        return false;
    } else if (distSquared > maxsq) {
        return false;
    } else {
        return true;
    }
}

function findConnections(points: AtomSpec[], otherPoints: AtomSpec[]) {
    for (let i = 0; i < points.length; i++) {
        const atom1 = points[i];
        for (let j = 0; j < otherPoints.length; j++) {
            const atom2 = otherPoints[j];

            if (areConnected(atom1, atom2)) {
                //gracefully handle one-sided bonds
                const a2i = atom1.bonds!.indexOf(atom2.index!);
                const a1i = atom2.bonds!.indexOf(atom1.index!);
                if (a2i === -1 && a1i === -1) {
                    atom1.bonds!.push(atom2.index!);
                    atom1.bondOrder!.push(1);
                    atom2.bonds!.push(atom1.index!);
                    atom2.bondOrder!.push(1);
                } else if (a2i === -1) {
                    atom1.bonds!.push(atom2.index!);
                    atom1.bondOrder!.push(atom2.bondOrder![a1i]);
                } else if (a1i === -1) {
                    atom2.bonds!.push(atom1.index!);
                    atom2.bondOrder!.push(atom1.bondOrder![a2i]);
                }
            }
        }
    }
}

export function assignBonds(atoms: AtomSpec[]): void {
    // assign bonds - yuck, can't count on connect records

    for (let i = 0, n = atoms.length; i < n; i++) {
        // Don't reindex if atoms are  already indexed
        if (!atoms[i].index) atoms[i].index = i;

        atoms[i].bonds = [];
        atoms[i].bondOrder = [];
    }

    const grid: Record<number, Record<number, Record<number, AtomSpec[]>>> = {};
    const MAX_BOND_LENGTH = 4.95; // (largest bond length, Cs) 2.25 * 2 * 1.1 (fudge factor)

    for (let index = 0; index < atoms.length; index++) {
        const atom = atoms[index];
        const x = Math.floor(atom.x! / MAX_BOND_LENGTH);
        const y = Math.floor(atom.y! / MAX_BOND_LENGTH);
        const z = Math.floor(atom.z! / MAX_BOND_LENGTH);
        if (!grid[x]) {
            grid[x] = {};
        }
        if (!grid[x][y]) {
            grid[x][y] = {};
        }
        if (!grid[x][y][z]) {
            grid[x][y][z] = [];
        }

        grid[x][y][z].push(atom);
    }

    const OFFSETS = [
        { x: 0, y: 0, z: 1 },
        { x: 0, y: 1, z: -1 },
        { x: 0, y: 1, z: 0 },
        { x: 0, y: 1, z: 1 },
        { x: 1, y: -1, z: -1 },
        { x: 1, y: -1, z: 0 },
        { x: 1, y: -1, z: 1 },
        { x: 1, y: 0, z: -1 },
        { x: 1, y: 0, z: 0 },
        { x: 1, y: 0, z: 1 },
        { x: 1, y: 1, z: -1 },
        { x: 1, y: 1, z: 0 },
        { x: 1, y: 1, z: 1 },
    ];

    for (const x_s in grid) {
        const x = parseInt(x_s, 10);
        for (const y_s in grid[x]) {
            const y = parseInt(y_s, 10);
            for (const z_s in grid[x][y]) {
                const z = parseInt(z_s, 10);
                const points = grid[x][y][z];

                for (let i = 0; i < points.length; i++) {
                    const atom1 = points[i];
                    for (let j = i + 1; j < points.length; j++) {
                        const atom2 = points[j];

                        if (areConnected(atom1, atom2)) {
                            if (atom1.bonds!.indexOf(atom2.index!) === -1) {
                                atom1.bonds!.push(atom2.index!);
                                atom1.bondOrder!.push(1);
                                atom2.bonds!.push(atom1.index!);
                                atom2.bondOrder!.push(1);
                            }
                        }
                    }
                }

                for (let o = 0; o < OFFSETS.length; o++) {
                    const offset = OFFSETS[o];
                    if (
                        !grid[x + offset.x] ||
                        !grid[x + offset.x][y + offset.y] ||
                        !grid[x + offset.x][y + offset.y][z + offset.z]
                    )
                        continue;

                    const otherPoints = grid[x + offset.x][y + offset.y][z + offset.z];
                    findConnections(points, otherPoints);
                }
            }
        }
    }
}

function assignBackboneHBonds(atomsarray: Array<AtomSpec>, hbondCutoff: number) {
    const maxlength = hbondCutoff || 3.2;
    const maxlengthSq = maxlength * maxlength;
    const atoms = [];

    for (let i = 0, n = atomsarray.length; i < n; i++) {
        atomsarray[i].index = i;
        // only consider 'N' and 'O'
        const atom = atomsarray[i];
        if (!atom.hetflag && (atom.atom === 'N' || atom.atom === 'O')) {
            atoms.push(atom);
            atom.hbondOther = null;
            atom.hbondDistanceSq = Number.POSITIVE_INFINITY;
        }
    }

    atoms.sort(function (a, b) {
        return a.z! - b.z!;
    });
    for (let i = 0, n = atoms.length; i < n; i++) {
        const ai = atoms[i];

        for (let j = i + 1; j < n; j++) {
            const aj = atoms[j];
            const zdiff = aj.z! - ai.z!;
            if (zdiff > maxlength)
                // can't be connected
                break;
            if (aj.atom === ai.atom) continue; // can't be connected, but later might be
            const ydiff = Math.abs(aj.y! - ai.y!);
            if (ydiff > maxlength) continue;
            const xdiff = Math.abs(aj.x! - ai.x!);
            if (xdiff > maxlength) continue;
            const dist = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;
            if (dist > maxlengthSq) continue;

            if (aj.chain === ai.chain && Math.abs(aj.resi! - ai.resi!) < 4) continue; // ignore bonds between too close residues
            // select closest hbond
            if (dist < ai.hbondDistanceSq!) {
                ai.hbondOther = aj;
                ai.hbondDistanceSq = dist;
            }
            if (dist < aj.hbondDistanceSq!) {
                aj.hbondOther = ai;
                aj.hbondDistanceSq = dist;
            }
        }
    }
}

export function computeSecondaryStructure(atomsarray: Array<AtomSpec>, hbondCutoff: number) {
    assignBackboneHBonds(atomsarray, hbondCutoff);

    // compute, per residue, what the secondary structure is
    const chres: Record<string, Array<string>> = {}; // lookup by chain and resid
    let i: number, il: number, c: string | number, r: number; // i: used in for loop, il: length of atomsarray
    let atom: AtomSpec, val: string;

    //identify helices first
    for (i = 0, il = atomsarray.length; i < il; i++) {
        atom = atomsarray[i];

        if (chres[atom.chain!] === undefined) {
            chres[atom.chain!] = [];
        }

        if (isFinite(atom.hbondDistanceSq!)) {
            const other = atom.hbondOther;
            if (chres[other.chain] === undefined) chres[other.chain] = [];

            if (Math.abs(other.resi - atom.resi!) === 4) {
                // helix
                chres[atom.chain!][atom.resi!] = 'h';
            }
        }
    }

    // plug gaps in helices
    for (c in chres) {
        for (r = 1; r < chres[c].length - 1; r++) {
            const valbefore = chres[c][r - 1];
            const valafter = chres[c][r + 1];
            val = chres[c][r];
            if (valbefore === 'h' && valbefore === valafter && val !== valbefore) {
                chres[c][r] = valbefore;
            }
        }
    }

    //now potential sheets - but only if mate not part of helix
    for (i = 0, il = atomsarray.length; i < il; i++) {
        atom = atomsarray[i];

        if (
            isFinite(atom.hbondDistanceSq!) &&
            chres[atom.chain!][atom.resi!] !== 'h' &&
            atom.ss !== 'h'
        ) {
            chres[atom.chain!][atom.resi!] = 'maybesheet';
        }
    }

    //sheets must bond to other sheets
    for (let i = 0, il = atomsarray.length; i < il; i++) {
        atom = atomsarray[i];

        if (isFinite(atom.hbondDistanceSq!) && chres[atom.chain!][atom.resi!] === 'maybesheet') {
            const other = atom.hbondOther;
            const otherval = chres[other.chain][other.resi];
            if (otherval === 'maybesheet' || otherval === 's') {
                // true sheet
                chres[atom.chain!][atom.resi!] = 's';
                chres[other.chain][other.resi] = 's';
            }
        }
    }

    // plug gaps in sheets and remove singletons
    for (const c in chres) {
        for (let r = 1; r < chres[c].length - 1; r++) {
            const valbefore = chres[c][r - 1];
            const valafter = chres[c][r + 1];
            val = chres[c][r];
            if (valbefore === 's' && valbefore === valafter && val !== valbefore) {
                chres[c][r] = valbefore;
            }
        }
        for (let r = 0; r < chres[c].length; r++) {
            const val = chres[c][r];
            if (val === 'h' || val === 's') {
                if (chres[c][r - 1] !== val && chres[c][r + 1] !== val) {
                    delete chres[c][r];
                }
            }
        }
    }

    // assign to all atoms in residue, keep track of start
    for (i = 0, il = atomsarray.length; i < il; i++) {
        atom = atomsarray[i];
        val = chres[atom.chain!][atom.resi!];

        //clear hbondOther to eliminate circular references that prohibit serialization
        delete atom.hbondOther;
        delete atom.hbondDistanceSq;
        if (val === undefined || val === 'maybesheet') continue;

        atom.ss = val;
        if (chres[atom.chain!][atom.resi! - 1] !== val) atom.ssbegin = true;
        if (chres[atom.chain!][atom.resi! + 1] !== val) atom.ssend = true;
    }
}
