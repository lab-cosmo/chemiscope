/**
 * @packageDocumentation
 * @module structure
 */

/// Cartoon decorations that 3Dmol does not draw by itself: nucleobases as
/// filled rings, and carbohydrates as 3D-SNFG symbols
/// (https://www.ncbi.nlm.nih.gov/glycans/snfg.html).

import * as $3Dmol from '3dmol';

type Vec = [number, number, number];

const add = (a: Vec, b: Vec): Vec => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const sub = (a: Vec, b: Vec): Vec => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const scale = (a: Vec, s: number): Vec => [a[0] * s, a[1] * s, a[2] * s];
const dot = (a: Vec, b: Vec): number => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a: Vec, b: Vec): Vec => [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
];
const normalize = (a: Vec): Vec => scale(a, 1 / (Math.sqrt(dot(a, a)) || 1));
const mean = (points: Vec[]): Vec => scale(points.reduce(add, [0, 0, 0]), 1 / points.length);
const position = (atom: $3Dmol.AtomSpec): Vec => [atom.x ?? 0, atom.y ?? 0, atom.z ?? 0];
const xyz = (v: Vec): $3Dmol.XYZ => ({ x: v[0], y: v[1], z: v[2] });

// ==================== residue classification ====================

const CHARMM_BASES: Record<string, string> = { ADE: 'A', CYT: 'C', GUA: 'G', THY: 'T', URA: 'U' };

/**
 * Get the nucleobase (A, C, G, T, U or I) of a nucleotide residue name, using
 * PDB (`DA`, `A`), Amber (`DA5`, `RA3`, `A5`) or CHARMM (`ADE`) conventions.
 */
export function nucleobase(resname: string): string | undefined {
    const name = resname.trim().toUpperCase();
    return CHARMM_BASES[name] ?? /^[DR]?([ACGTUI])[35N]?$/.exec(name)?.[1];
}

/**
 * Residue name to give to 3Dmol, which only draws nucleic acid cartoons for
 * the PDB names of the standard nucleotides.
 */
export function cartoonResname(resname: string): string {
    const base = nucleobase(resname);
    return base === undefined ? resname : base === 'I' ? 'G' : base;
}

/** 2D shapes of the SNFG, drawn as the 3D shapes of 3D-SNFG */
type SnfgShape =
    | 'filled-circle'
    | 'filled-square'
    | 'crossed-square'
    | 'divided-diamond'
    | 'filled-triangle'
    | 'divided-triangle'
    | 'flat-rectangle'
    | 'filled-star'
    | 'filled-diamond'
    | 'flat-diamond'
    | 'flat-hexagon'
    | 'pentagon';

/** RGB colors from table 2 of https://www.ncbi.nlm.nih.gov/glycans/snfg.html */
const SNFG_COLORS = {
    white: 0xffffff,
    blue: 0x0072bc,
    green: 0x00a651,
    yellow: 0xffd400,
    orange: 0xf47920,
    pink: 0xf69ea1,
    purple: 0xa54399,
    lightBlue: 0x8fcce9,
    brown: 0xa17a4d,
    red: 0xed1c24,
};

/**
 * Table 1 of https://www.ncbi.nlm.nih.gov/glycans/snfg.html, without the
 * generic (white) column. `-` marks empty cells.
 */
// prettier-ignore
const SNFG_TABLE: [SnfgShape, string][] = [
    //                   blue    green    yellow  orange  pink     purple  light blue brown red
    ['filled-circle',   'Glc     Man      Gal     Gul     Alt      All     Tal      Ido     -'],
    ['filled-square',   'GlcNAc  ManNAc   GalNAc  GulNAc  AltNAc   AllNAc  TalNAc   IdoNAc  -'],
    ['crossed-square',  'GlcN    ManN     GalN    GulN    AltN     AllN    TalN     IdoN    -'],
    ['divided-diamond', 'GlcA    ManA     GalA    GulA    AltA     AllA    TalA     IdoA    -'],
    ['filled-triangle', 'Qui     Rha      -       6dGul   6dAlt    -       6dTal    -       Fuc'],
    ['divided-triangle','QuiNAc  RhaNAc   -       -       6dAltNAc -       6dTalNAc -       FucNAc'],
    ['flat-rectangle',  'Oli     Tyv      -       Abe     Par      Dig     Col      -       -'],
    ['filled-star',     '-       Ara      Lyx     Xyl     Rib      -       -        -       -'],
    ['filled-diamond',  '-       Kdn      -       -       -        Neu5Ac  Neu5Gc   Neu     Sia'],
    ['flat-diamond',    '-       Pse      Leg     -       Aci      -       4eLeg    -       -'],
    ['flat-hexagon',    'Bac     LDmanHep Kdo     Dha     DDmanHep MurNAc  MurNGc   Mur     -'],
    ['pentagon',        'Api     Fru      Tag     Sor     Psi      -       -        -       -'],
];

/** SNFG shape and color for each monosaccharide */
const SNFG = new Map<string, [SnfgShape, number]>();
for (const [shape, row] of SNFG_TABLE) {
    const { blue, green, yellow, orange, pink, purple, lightBlue, brown, red } = SNFG_COLORS;
    const colors = [blue, green, yellow, orange, pink, purple, lightBlue, brown, red];
    row.split(/\s+/).forEach((name, i) => name !== '-' && SNFG.set(name, [shape, colors[i]]));
}

/** L-hexuronic acids have the colored half of the diamond at the bottom */
const BOTTOM_COLORED = new Set(['AltA', 'IdoA']);

/** PDB chemical components and CHARMM residue names of monosaccharides */
const SACCHARIDE_NAMES: Record<string, string> = {};
for (const [saccharide, names] of Object.entries({
    Glc: 'GLC BGC AGLC BGLC',
    Man: 'MAN BMA AMAN BMAN',
    Gal: 'GLA GAL AGAL BGAL',
    GlcNAc: 'NAG NDG AGLCNA BGLCNA BGLCN0',
    ManNAc: 'BM3 BM7',
    GalNAc: 'A2G NGA AGALNA BGALNA',
    GlcN: 'PA1 GCS',
    GalN: 'X6X 1GN',
    GlcA: 'GCU BDP AGLCA BGLCA BGLCA0',
    ManA: 'MAV BEM',
    GalA: 'ADA GTR',
    IdoA: 'IDR AIDOA BIDOA',
    Qui: 'G6D',
    Rha: 'RAM RM4 ARHM BRHM',
    Fuc: 'FUC FUL AFUC BFUC',
    FucNAc: '49T',
    Tyv: 'TYV',
    Abe: 'ABE',
    Ara: 'ARA ARB AHR FUB AARB BARB',
    Xyl: 'XYS XYP AXYL BXYL',
    Rib: 'RIB ARIB BRIB',
    Kdn: 'KDN KDM',
    Neu5Ac: 'SIA SLB ANE5AC BNE5AC',
    Neu5Gc: 'NGC NGE',
    Kdo: 'KDO',
    LDmanHep: 'GMH',
    MurNAc: 'MUB AMU',
    Fru: 'FRU BDF AFRU BFRU',
})) {
    names.split(' ').forEach((name) => (SACCHARIDE_NAMES[name] = saccharide));
}

/** Second letter of GLYCAM residue names, which are lower case for L sugars */
// prettier-ignore
const GLYCAM_LETTERS: Record<string, string> = {
    A: 'Ara', B: 'Sor', C: 'Fru', D: 'Lyx', E: 'Alt', F: 'Fuc', G: 'Glc', H: 'Rha', I: 'Ido',
    J: 'Tag', K: 'Gul', L: 'Gal', M: 'Man', N: 'All', O: 'GalA', P: 'Psi', Q: 'Qui', R: 'Rib',
    S: 'Neu5Ac', T: 'Tal', U: 'IdoA', V: 'GalNAc', W: 'ManNAc', X: 'Xyl', Y: 'GlcNAc', Z: 'GlcA',
};

/**
 * Get the monosaccharide (e.g. `GlcNAc`) of a residue name, using PDB, CHARMM
 * or GLYCAM conventions. GLYCAM names (`<linkage><sugar><anomer>`, e.g. `4YB`)
 * can collide with other residue names, so the ring should be checked with
 * {@link saccharideRing} as well.
 */
export function saccharide(resname: string): string | undefined {
    const name = resname.trim();
    const glycam = /^[0-9A-Z]([A-Za-z])[AB]$/.exec(name);
    return (
        SACCHARIDE_NAMES[name.toUpperCase()] ??
        (glycam === null ? undefined : GLYCAM_LETTERS[glycam[1].toUpperCase()])
    );
}

/** Candidate rings, starting with the anomeric carbon */
const SACCHARIDE_RINGS = [
    ['C1', 'C2', 'C3', 'C4', 'C5', 'O5'], // aldopyranoses
    ['C2', 'C3', 'C4', 'C5', 'C6', 'O6'], // ketopyranoses, Neu5Ac, Kdo
    ['C1', 'C2', 'C3', 'C4', 'O4'], // aldofuranoses
    ['C2', 'C3', 'C4', 'C5', 'O5'], // ketofuranoses
];

/**
 * Find the ring of a carbohydrate residue, as a list of atoms starting with the
 * anomeric carbon. Returns `undefined` if this is not a known monosaccharide.
 */
export function saccharideRing(residue: $3Dmol.AtomSpec[]): $3Dmol.AtomSpec[] | undefined {
    if (residue.length === 0 || saccharide(residue[0].resn ?? '') === undefined) {
        return undefined;
    }
    const byName = new Map(residue.map((atom) => [atom.atom, atom]));
    for (const names of SACCHARIDE_RINGS) {
        const ring = names.map((name) => byName.get(name));
        if (ring.every((atom) => atom !== undefined) && isRing(ring.map(position))) {
            return ring;
        }
    }
    return undefined;
}

function isRing(points: Vec[]): boolean {
    return points.every((p, i) => {
        const d = sub(p, points[(i + 1) % points.length]);
        return dot(d, d) < 1.9 * 1.9;
    });
}

/** Group consecutive atoms belonging to the same residue */
export function residues(atoms: $3Dmol.AtomSpec[]): $3Dmol.AtomSpec[][] {
    const groups: $3Dmol.AtomSpec[][] = [];
    let previous: $3Dmol.AtomSpec | undefined;
    for (const atom of atoms) {
        if (
            previous === undefined ||
            atom.chain !== previous.chain ||
            atom.resi !== previous.resi ||
            atom.resn !== previous.resn
        ) {
            groups.push([]);
        }
        groups[groups.length - 1].push(atom);
        previous = atom;
    }
    return groups;
}

// ==================== meshes ====================

/** 3Dmol custom shapes can not hold more vertices than this */
const MAX_VERTICES = 60000;

interface MeshChunk {
    vertexArr: $3Dmol.XYZ[];
    normalArr: $3Dmol.XYZ[];
    faceArr: number[];
    color: number[];
}

/**
 * Triangle mesh with per-vertex colors, split in chunks small enough for
 * 3Dmol custom shapes
 */
class Mesh {
    public readonly chunks: MeshChunk[] = [];
    private _current: MeshChunk = this._newChunk();

    private _newChunk(): MeshChunk {
        const chunk = { vertexArr: [], normalArr: [], faceArr: [], color: [] };
        this.chunks.push(chunk);
        return chunk;
    }

    /** Make sure the current chunk has room for `count` more vertices */
    private _reserve(count: number): MeshChunk {
        if (this._current.vertexArr.length + count > MAX_VERTICES) {
            this._current = this._newChunk();
        }
        return this._current;
    }

    private _vertex(point: Vec, normal: Vec, color: number): number {
        const chunk = this._current;
        chunk.vertexArr.push(xyz(point));
        chunk.normalArr.push(xyz(normal));
        chunk.color.push(color);
        return chunk.vertexArr.length - 1;
    }

    /** Add a flat triangle, facing away from the `inside` point */
    public triangle(a: Vec, b: Vec, c: Vec, inside: Vec, color: number): void {
        let normal = normalize(cross(sub(b, a), sub(c, a)));
        if (dot(normal, sub(a, inside)) < 0) {
            [b, c] = [c, b];
            normal = scale(normal, -1);
        }
        this._reserve(3).faceArr.push(
            this._vertex(a, normal, color),
            this._vertex(b, normal, color),
            this._vertex(c, normal, color)
        );
    }

    /** Add a convex polygon, as a fan around its center */
    public polygon(points: Vec[], inside: Vec, color: number): void {
        const center = mean(points);
        points.forEach((p, i) =>
            this.triangle(center, p, points[(i + 1) % points.length], inside, color)
        );
    }

    /** Extrude the polygon `points` by `±offset`, it must be star-shaped around its center */
    public slab(points: Vec[], offset: Vec, color: number): void {
        const center = mean(points);
        const top = points.map((p) => add(p, offset));
        const bottom = points.map((p) => sub(p, offset));
        this.polygon(top, center, color);
        this.polygon(bottom, center, color);
        points.forEach((_, i) => {
            const j = (i + 1) % points.length;
            this.triangle(top[i], top[j], bottom[j], center, color);
            this.triangle(top[i], bottom[j], bottom[i], center, color);
        });
    }

    /** Pyramid with the given (convex) base and apex */
    public pyramid(base: Vec[], apex: Vec, color: number): void {
        const inside = mean([...base, apex]);
        this.polygon(base, inside, color);
        base.forEach((p, i) => this.triangle(p, base[(i + 1) % base.length], apex, inside, color));
    }

    /** Smooth surface from a grid of `rows x columns` points, periodic along columns */
    private _surface(
        rows: number,
        columns: number,
        point: (i: number, j: number) => [Vec, Vec],
        color: number
    ): void {
        const chunk = this._reserve(rows * columns);
        const start = chunk.vertexArr.length;
        for (let i = 0; i < rows; i++) {
            for (let j = 0; j < columns; j++) {
                this._vertex(...point(i, j), color);
            }
        }
        for (let i = 0; i + 1 < rows; i++) {
            for (let j = 0; j < columns; j++) {
                const a = start + i * columns + j;
                const b = start + i * columns + ((j + 1) % columns);
                chunk.faceArr.push(a, a + columns, b, b, a + columns, b + columns);
            }
        }
    }

    public sphere(center: Vec, radius: number, color: number): void {
        this._surface(
            11,
            16,
            (i, j) => {
                const [theta, phi] = [(Math.PI * i) / 10, (2 * Math.PI * j) / 16];
                const n: Vec = [
                    Math.sin(theta) * Math.cos(phi),
                    Math.sin(theta) * Math.sin(phi),
                    Math.cos(theta),
                ];
                return [add(center, scale(n, radius)), n];
            },
            color
        );
    }

    /** Cylinder with flat caps */
    public cylinder(start: Vec, end: Vec, radius: number, color: number): void {
        const axis = normalize(sub(end, start));
        const p = normalize(cross(axis, Math.abs(axis[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]));
        const q = cross(axis, p);
        const columns = 12;
        const radial = (j: number): Vec => {
            const angle = (2 * Math.PI * j) / columns;
            return add(scale(p, Math.cos(angle)), scale(q, Math.sin(angle)));
        };
        this._surface(
            2,
            columns,
            (i, j) => [add(i === 0 ? end : start, scale(radial(j), radius)), radial(j)],
            color
        );
        const middle = mean([start, end]);
        for (const center of [start, end]) {
            const cap = Array.from({ length: columns }, (_, j) =>
                add(center, scale(radial(j), radius))
            );
            this.polygon(cap, middle, color);
        }
    }
}

/** Orthonormal frame centered on a saccharide ring */
interface Frame {
    center: Vec;
    /** in the ring plane, towards the anomeric carbon */
    u: Vec;
    /** in the ring plane, orthogonal to `u` */
    v: Vec;
    /** normal to the ring plane */
    n: Vec;
}

function ringNormal(points: Vec[]): Vec {
    const center = mean(points);
    return normalize(
        points.reduce<Vec>(
            (n, p, i) =>
                add(n, cross(sub(p, center), sub(points[(i + 1) % points.length], center))),
            [0, 0, 0]
        )
    );
}

function ringFrame(points: Vec[]): Frame {
    const center = mean(points);
    const n = ringNormal(points);
    const toAnomeric = sub(points[0], center);
    const u = normalize(sub(toAnomeric, scale(n, dot(toAnomeric, n))));
    return { center, u, v: cross(n, u), n };
}

/** Size of the SNFG symbols, in Angstrom */
const SNFG_SIZE = 1.6;

/**
 * Add the 3D version of a SNFG symbol to the mesh. The 2D symbol lies in the
 * (u, v) plane of the ring with v pointing up, using the proportions of the
 * symbols on the SNFG page.
 */
function addSnfgSymbol(
    mesh: Mesh,
    frame: Frame,
    shape: SnfgShape,
    color: number,
    bottomColored: boolean
): void {
    const { center, u, v, n } = frame;
    const { white } = SNFG_COLORS;
    const R = SNFG_SIZE;
    // local (x, y, z) coordinates along (u, v, n) to world coordinates
    const at = (x: number, y: number, z = 0): Vec =>
        add(center, add(scale(u, x), add(scale(v, y), scale(n, z))));
    // regular polygon/star with a vertex pointing up
    const polygon = (count: number, radius: (i: number) => number) =>
        Array.from({ length: count }, (_, i) => {
            const angle = Math.PI / 2 + (2 * Math.PI * i) / count;
            return at(radius(i) * Math.cos(angle), radius(i) * Math.sin(angle));
        });
    const flat = (points: Vec[], thickness: number, c = color) =>
        mesh.slab(points, scale(n, thickness / 2), c);
    // octahedron, split into top (y > 0) and bottom (y < 0) halves
    const diamond = (topColor: number, bottomColor: number) => {
        const a = 1.25 * R;
        const square = [at(a, 0), at(0, 0, a), at(-a, 0), at(0, 0, -a)];
        mesh.pyramid(square, at(0, a), topColor);
        mesh.pyramid(square, at(0, -a), bottomColor);
    };
    // cone pointing up, split into right (x > 0) and left (x < 0) halves
    const cone = (rightColor: number, leftColor: number) => {
        const h = 1.9 * R;
        const half = (sign: number) =>
            Array.from({ length: 9 }, (_, i) => {
                const angle = (Math.PI * i) / 8 - Math.PI / 2;
                return at((sign * h * Math.cos(angle)) / 2, -h / 2, (h * Math.sin(angle)) / 2);
            });
        mesh.pyramid(half(1), at(0, h / 2), rightColor);
        mesh.pyramid(half(-1), at(0, h / 2), leftColor);
    };

    const e = 0.8 * R; // half edge of the cubes
    switch (shape) {
        case 'filled-circle':
            mesh.sphere(center, R, color);
            break;
        case 'filled-square':
            flat([at(-e, -e), at(e, -e), at(e, e), at(-e, e)], 2 * e);
            break;
        case 'crossed-square':
            // colored upper right half, white lower left half
            flat([at(-e, e), at(e, e), at(e, -e)], 2 * e);
            flat([at(-e, e), at(-e, -e), at(e, -e)], 2 * e, white);
            break;
        case 'divided-diamond':
            diamond(bottomColored ? white : color, bottomColored ? color : white);
            break;
        case 'filled-diamond':
            diamond(color, color);
            break;
        case 'filled-triangle':
            cone(color, color);
            break;
        case 'divided-triangle':
            // colored right half, white left half
            cone(color, white);
            break;
        case 'flat-rectangle': {
            const [w, h] = [1.25 * R, (1.25 * R) / Math.sqrt(3)];
            flat([at(-w, -h), at(w, -h), at(w, h), at(-w, h)], 2 * h);
            break;
        }
        case 'filled-star':
            flat(
                polygon(10, (i) => (i % 2 === 0 ? 1.35 : 0.382 * 1.35) * R),
                0.6 * R
            );
            break;
        case 'flat-diamond': {
            const [w, h] = [1.35 * R, 0.664 * 1.35 * R];
            flat([at(w, 0), at(0, h), at(-w, 0), at(0, -h)], 0.6 * R);
            break;
        }
        case 'flat-hexagon': {
            const [w, h] = [1.3 * R, 0.612 * 1.3 * R];
            flat(
                [at(w, 0), at(w / 2, h), at(-w / 2, h), at(-w, 0), at(-w / 2, -h), at(w / 2, -h)],
                0.6 * R
            );
            break;
        }
        case 'pentagon':
            flat(
                polygon(5, () => 1.25 * R),
                0.6 * R
            );
            break;
    }
}

// ==================== decorations ====================

/** Base colors following the NDB/DSSR convention */
const BASE_COLORS: Record<string, number> = {
    A: 0xe03c31,
    C: 0xf5c400,
    G: 0x2ca02c,
    T: 0x1f5fc8,
    U: 0x17becf,
};

const PYRIMIDINE_RING = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'];
const IMIDAZOLE_RING = ['C4', 'C5', 'N7', 'C8', 'N9'];
/** Thickness of the base slabs and radius of the links, in Angstrom */
const [BASE_THICKNESS, LINK_RADIUS, LINK_COLOR] = [0.5, 0.3, 0xaaaaaa];

export interface Decorations {
    /** meshes for base rings, saccharide symbols and the links between them */
    meshes: $3Dmol.CustomShapeSpec[];
    /** indexes of atoms drawn as a saccharide symbol instead */
    hidden: Set<number>;
}

/**
 * Compute nucleobase rings (if `cartoon` is true) and SNFG symbols for
 * carbohydrates (if `snfg` is true) for the given 3Dmol `atoms`, which should
 * already have bonds and residue information.
 */
export function biomolDecorations(
    atoms: $3Dmol.AtomSpec[],
    cartoon: boolean,
    snfg: boolean
): Decorations {
    const mesh = new Mesh();
    const hidden = new Set<number>();
    const groups = residues(atoms);

    if (cartoon) {
        groups.forEach((residue, i) => {
            const base = nucleobase(residue[0].resn ?? '');
            if (base === undefined || residue[0].hetflag) {
                return;
            }
            const byName = new Map(residue.map((atom) => [atom.atom ?? '', atom]));
            const ring = (names: string[]) => names.map((name) => byName.get(name));
            const pyrimidine = ring(PYRIMIDINE_RING);
            const imidazole = ring(IMIDAZOLE_RING);
            const anchor = backboneAnchor(residue, groups[i + 1]);
            if (pyrimidine.some((atom) => atom === undefined) || anchor === undefined) {
                return;
            }
            const color = BASE_COLORS[base] ?? BASE_COLORS.G;
            const purine = imidazole.every((atom) => atom !== undefined);
            const rings = purine ? [pyrimidine, imidazole] : [pyrimidine];
            for (const atomsInRing of rings) {
                const points = atomsInRing.map((atom) => position(atom as $3Dmol.AtomSpec));
                mesh.slab(points, scale(ringNormal(points), BASE_THICKNESS / 2), color);
            }
            // glycosidic nitrogen: N9 for purines, N1 for pyrimidines
            const glycosidic = (purine ? imidazole[4] : pyrimidine[0]) as $3Dmol.AtomSpec;
            mesh.cylinder(anchor, position(glycosidic), LINK_RADIUS, color);
        });
    }

    if (snfg) {
        const residueOf = new Map<number | undefined, number>();
        groups.forEach((residue, i) => residue.forEach((atom) => residueOf.set(atom.index, i)));
        const byIndex = new Map(atoms.map((atom) => [atom.index, atom]));
        const heavy = atoms.filter((atom) => atom.elem !== 'H');

        const symbols = new Map<number, { center: Vec; anomeric: $3Dmol.AtomSpec }>();
        groups.forEach((residue, i) => {
            const ring = saccharideRing(residue);
            if (ring === undefined) {
                return;
            }
            const name = saccharide(residue[0].resn ?? '') ?? '';
            // unknown monosaccharides are white flat hexagons
            const [shape, color] = SNFG.get(name) ?? ['flat-hexagon', SNFG_COLORS.white];
            const frame = ringFrame(ring.map(position));
            addSnfgSymbol(mesh, frame, shape, color, BOTTOM_COLORED.has(name));
            symbols.set(i, { center: frame.center, anomeric: ring[0] });
            residue.forEach((atom) => hidden.add(atom.index ?? -1));
        });

        const linked = new Set<string>();
        for (const [i, { center, anomeric }] of symbols) {
            const outside = (atom: $3Dmol.AtomSpec | undefined): atom is $3Dmol.AtomSpec =>
                atom !== undefined && residueOf.get(atom.index) !== i;
            let partners = groups[i]
                .flatMap((atom) => (atom.bonds ?? []).map((j) => byIndex.get(j)))
                .filter(outside);
            if (partners.length === 0) {
                // no bonds between residues in the input, look for the
                // glycosidic bond of the anomeric carbon instead
                const p = position(anomeric);
                partners = heavy.filter((atom) => {
                    const d = sub(position(atom), p);
                    return outside(atom) && dot(d, d) < 1.9 * 1.9;
                });
            }

            for (const partner of partners) {
                const j = residueOf.get(partner.index) as number;
                const key = i < j ? `${i}-${j}` : `${j}-${i}`;
                if (linked.has(key)) {
                    continue;
                }
                linked.add(key);
                // link to the other symbol, or to the alpha carbon of
                // proteins drawn as cartoon, or else to the bonded atom
                const ca =
                    cartoon && !partner.hetflag
                        ? groups[j].find((atom) => atom.atom === 'CA' && atom.elem === 'C')
                        : undefined;
                const end = symbols.get(j)?.center ?? position(ca ?? partner);
                mesh.cylinder(center, end, LINK_RADIUS, LINK_COLOR);
            }
        }
    }

    return { meshes: mesh.chunks.filter((chunk) => chunk.faceArr.length > 0), hidden };
}

/**
 * Point on the nucleic acid backbone ribbon to which the base is attached,
 * using the same convention as 3Dmol: the middle of the phosphorus atoms of
 * this residue and the next one.
 */
function backboneAnchor(
    residue: $3Dmol.AtomSpec[],
    next: $3Dmol.AtomSpec[] | undefined
): Vec | undefined {
    const find = (atoms: $3Dmol.AtomSpec[] | undefined, names: string[]) =>
        atoms?.find((atom) => names.includes(atom.atom ?? ''));
    const backbone = find(residue, ['P']) ?? find(residue, ["O5'", 'O5*']);
    if (backbone === undefined) {
        return undefined;
    }
    const nextP = next?.[0].chain === residue[0].chain ? find(next, ['P']) : undefined;
    const end = nextP ?? find(residue, ["O3'", 'O3*']);
    if (end === undefined) {
        return position(backbone);
    }
    const [a, b] = [position(backbone), position(end)];
    const d = sub(b, a);
    // 3Dmol breaks the ribbon for phosphorus atoms further apart than 8 A
    return dot(d, d) < 64 ? mean([a, b]) : a;
}
