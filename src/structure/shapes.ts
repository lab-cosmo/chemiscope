/**
 *
 * This file consists of shapes supported by chemiscope via 3dmol
 *
 * To add a new shape type, extend the `Shape` class and define the inherited
 * functions (`validateParams` and `outputTo3Dmol`). `validateParams` should
 * check a `ShapeData` object to see that it has the correct parameters for the
 * given type. `outputTo3Dmol` should take a color and output a
 * `CustomShapeSpec` (https://3dmol.csb.pitt.edu/doc/types.html#CustomShapeSpec)
 * object. You will also need to add any parameters for your shape to ShapeData.
 *
 * To add a new shape, you will also need to add code to `structure/viewer.ts`.
 */
import assert from 'assert';

import { default as $3Dmol } from '3dmol';
import { Quaternion } from '3dmol';
type ColorSpec = number | string | { r: number; g: number; b: number; a?: number }; // mysterious compile errors if I import from 3dmol

// Just an interface to enforce XYZ-type coordinates
export interface XYZ {
    x: number;
    y: number;
    z: number;
}

export interface BaseShapeData {
    position: [number, number, number];
    orientation?: [number, number, number, number];
    scale?: number;
    color?: ColorSpec;
}

export interface BaseShapeParameters<T> {
    kind: string;
    parameters: {
        global: Partial<T>;
        structure?: Array<Partial<T>>;
        atom?: Array<Partial<T>>;
    };
}

export interface SphereData extends BaseShapeData {
    radius: number;
}

export interface SphereParameters extends BaseShapeParameters<SphereData> {
    kind: 'sphere';
}

// Interface for ellipsoidal data, where the shape
// of the ellipsoid is given by the semiaxes, and the orientation
// by the general `orientation` property
export interface EllipsoidData extends BaseShapeData {
    semiaxes: [number, number, number];
}

export interface EllipsoidParameters extends BaseShapeParameters<EllipsoidData> {
    kind: 'ellipsoid';
}

// Interface for arrow data
export interface ArrowData extends BaseShapeData {
    vector: [number, number, number];
    base_radius?: number;
    head_radius?: number;
    head_length?: number;
}

export interface ArrowParameters extends BaseShapeParameters<ArrowData> {
    kind: 'arrow';
}

// Interface for polytope data, where
// orientation is stored in the (w, x, y, z) convention
// and simplices refers to the indices of the facets
export interface CustomShapeData extends BaseShapeData {
    vertices: [number, number, number][];
    simplices: [number, number, number][];
}

export interface CustomShapeParameters extends BaseShapeParameters<CustomShapeData> {
    kind: 'custom';
}

export type ShapeData = SphereData | EllipsoidData | ArrowData | CustomShapeData;
export type ShapeParameters =
    | SphereParameters
    | EllipsoidParameters
    | ArrowParameters
    | CustomShapeParameters;

function addXYZ(a: XYZ, b: XYZ): XYZ {
    return { x: a.x + b.x, y: a.y + b.y, z: a.z + b.z };
}

function subXYZ(a: XYZ, b: XYZ): XYZ {
    return { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z };
}

function multXYZ(a: XYZ, b: number): XYZ {
    return { x: a.x * b, y: a.y * b, z: a.z * b };
}

function dotXYZ(a: XYZ, b: XYZ): number {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

function crossXYZ(a: XYZ, b: XYZ): XYZ {
    return {
        x: a.y * b.z - a.z * b.y,
        y: a.z * b.x - a.x * b.z,
        z: a.x * b.y - a.y * b.x,
    };
}

function rotateAndPlace(vertex: XYZ, quaternion: Quaternion, position: XYZ): XYZ {
    const vertexQuaternion = new Quaternion(vertex.x, vertex.y, vertex.z, 0);
    vertexQuaternion.multiplyQuaternions(quaternion.clone(), vertexQuaternion);
    vertexQuaternion.multiply(quaternion.clone().inverse());

    return addXYZ(
        { x: vertexQuaternion.x, y: vertexQuaternion.y, z: vertexQuaternion.z },
        position
    );
}

function determineNormals(vertices: XYZ[], simplices: [number, number, number][]): XYZ[] {
    const vertexNormals: XYZ[] = [];
    const nFaces: number[] = [];

    for (let i = 0; i < vertices.length; i++) {
        vertexNormals.push({ x: 0, y: 0, z: 0 });
        nFaces.push(0);
    }

    for (const s of simplices) {
        const faceNormal: XYZ = crossXYZ(
            subXYZ(vertices[s[1]], vertices[s[0]]),
            subXYZ(vertices[s[2]], vertices[s[0]])
        );

        for (const ss of s) {
            vertexNormals[ss] = addXYZ(vertexNormals[ss], faceNormal);
            nFaces[ss] += 1;
        }
    }

    for (let i = 0; i < vertices.length; i++) {
        vertexNormals[i].x /= nFaces[i];
        vertexNormals[i].y /= nFaces[i];
        vertexNormals[i].z /= nFaces[i];
    }

    return vertexNormals;
}

export class Shape {
    public position: XYZ;
    /// orientation of the particle
    public orientation: Quaternion;
    /// scaling factor
    public scale: number;

    // orientation is passed to 3dmol in the (x, y, z, w) convention
    constructor(data: Partial<SphereData>) {
        const [qx, qy, qz, qw] = data.orientation || [0, 0, 0, 1];
        this.orientation = new Quaternion(qx, qy, qz, qw);
        const [x, y, z] = data.position || [0, 0, 0];
        this.position = { x, y, z };
        this.scale = data.scale || 1.0;
    }

    // disabling eslint because it complains parameters is never used
    // function to check the parameters corresponding to the given
    // sub-class exist.
    /* eslint-disable-next-line */
    public static validateParams(parameters: Record<string, unknown>): string {
        return '';
    }
    // function to output a `CustomShapeSpec` for the given object
    public outputTo3Dmol(color: $3Dmol.ColorSpec): $3Dmol.CustomShapeSpec {
        return { color: color };
    }
}

function triangulateEllipsoid(
    semiaxes: [number, number, number],
    resolution: number = 20
): { vertices: XYZ[]; indices: number[] } {
    const step_eta = 360 / resolution;
    const step_omega = 180 / resolution;
    const [a1, a2, a3] = semiaxes;

    const indices: number[] = [];
    const vertices: XYZ[] = [];

    let eta, omega;
    let eta_rad, omega_rad;
    let index = 0;

    for (eta = -90; eta <= 90; eta += step_eta) {
        eta_rad = (eta * Math.PI) / 180.0;

        for (omega = -180; omega <= 180; omega += step_omega, index += 3) {
            omega_rad = (omega * Math.PI) / 180.0;

            const newVertex: XYZ = {
                x: a1 * Math.cos(eta_rad) * Math.cos(omega_rad),
                y: a2 * Math.cos(eta_rad) * Math.sin(omega_rad),
                z: a3 * Math.sin(eta_rad),
            };
            vertices.push(newVertex);
        }
    }

    // Now triangulation. I believe that it needs a counterclockwise
    // specification of vertices (else I have back-face) each triangle is
    // specified with 3 indices (gl.TRIANGLES).
    const num_samples_omega = 360 / step_omega;
    const num_samples_eta = 180 / step_eta;

    let i, j;
    for (i = 0, index = 0; i < num_samples_eta; i++) {
        for (j = 0; j < num_samples_omega; j++, index += 6) {
            indices[index] = i * num_samples_omega + i + j;
            indices[index + 1] = (i + 1) * num_samples_omega + 2 + i + j;
            indices[index + 2] = (i + 1) * num_samples_omega + 1 + i + j;

            indices[index + 3] = i * num_samples_omega + i + j;
            indices[index + 4] = i * num_samples_omega + 1 + i + j;
            indices[index + 5] = (i + 1) * num_samples_omega + 2 + i + j;
        }
    }
    return {
        vertices: vertices,
        indices: indices,
    };
}

function validateOrientation(orientation: unknown, kind: string): string {
    if (!Array.isArray(orientation) || orientation.length !== 4) {
        return `"orientation" must be an array with 4 elements for "${kind}" shapes`;
    }

    const [qx, qy, qz, qw] = orientation as unknown[];
    if (
        typeof qx !== 'number' ||
        typeof qy !== 'number' ||
        typeof qz !== 'number' ||
        typeof qw !== 'number'
    ) {
        return `"orientation" elements must be numbers for "${kind}" shapes`;
    }

    const norm2 = qx * qx + qy * qy + qz * qz + qw * qw;
    if (Math.abs(norm2 - 1) > 1e-6) {
        throw Error(`orientation must be normalized to 1 for "${kind}" shapes`);
    }

    return '';
}

function isPositiveInteger(o: unknown): boolean {
    return typeof o === 'number' && Number.isInteger(o) && o >= 0;
}

export class Sphere extends Shape {
    public radius: number;

    constructor(data: Partial<SphereData>) {
        super(data);
        this.radius = (data.radius || 1.0) * this.scale;
    }

    public static validateParameters(parameters: Record<string, unknown>): string {
        if (!('radius' in parameters)) {
            return '"radius" is required for "sphere" shapes';
        }

        if (typeof parameters.radius !== 'number') {
            return '"radius" must be a number in "sphere" shape';
        }

        return '';
    }

    public outputTo3Dmol(color: $3Dmol.ColorSpec, resolution: number = 20): $3Dmol.CustomShapeSpec {
        // Following https://computergraphics.stackexchange.com/questions/7426/triangulation-of-vertices-of-an-ellipsoid
        const triangulation = triangulateEllipsoid(
            [this.radius, this.radius, this.radius],
            resolution
        );
        const rawVertices = triangulation.vertices;
        const indices = triangulation.indices;
        const vertices: XYZ[] = [];
        const normals: XYZ[] = [];

        for (const v of rawVertices) {
            const newVertex: XYZ = addXYZ(v, this.position);
            vertices.push(newVertex);

            const newNormal: XYZ = {
                x: v.x / Math.pow(this.radius, 2.0),
                y: v.y / Math.pow(this.radius, 2.0),
                z: v.z / Math.pow(this.radius, 2.0),
            };
            normals.push(newNormal);
        }

        return {
            vertexArr: vertices,
            normalArr: normals,
            faceArr: indices,
            color: color,
        };
    }
}

export class Ellipsoid extends Shape {
    public semiaxes: [number, number, number];

    constructor(data: Partial<EllipsoidData>) {
        super(data);
        assert(data.semiaxes);
        this.semiaxes = [
            this.scale * data.semiaxes[0],
            this.scale * data.semiaxes[1],
            this.scale * data.semiaxes[2],
        ];
    }

    public static validateParameters(parameters: Record<string, unknown>): string {
        if (!('semiaxes' in parameters)) {
            return '"semiaxes" is required for "ellipsoid" shapes';
        }

        if (!Array.isArray(parameters.semiaxes) || parameters.semiaxes.length !== 3) {
            return '"semiaxes" must be an array with 3 elements for "ellipsoid" shapes';
        }

        const [ax, ay, az] = parameters.semiaxes as unknown[];
        if (typeof ax !== 'number' || typeof ay !== 'number' || typeof az !== 'number') {
            return '"semiaxes" elements must be numbers for "ellipsoid" shapes';
        }

        if ('orientation' in parameters) {
            const message = validateOrientation(parameters.orientation, 'ellipsoid');
            if (message !== '') {
                return message;
            }
        }

        return '';
    }

    public outputTo3Dmol(color: $3Dmol.ColorSpec, resolution: number = 20): $3Dmol.CustomShapeSpec {
        const triangulation = triangulateEllipsoid(this.semiaxes, resolution);
        const rawVertices = triangulation.vertices;
        const indices = triangulation.indices;
        const vertices: XYZ[] = [];
        const normals: XYZ[] = [];

        for (const v of rawVertices) {
            const newVertex: XYZ = rotateAndPlace(v, this.orientation, this.position);
            vertices.push(newVertex);

            // the normal to the ellipsoid surface is best computed in the original
            // coordinate system, and then rotated in place
            const newNormal: XYZ = rotateAndPlace(
                {
                    x: v.x / Math.pow(this.semiaxes[0], 2.0),
                    y: v.y / Math.pow(this.semiaxes[1], 2.0),
                    z: v.z / Math.pow(this.semiaxes[2], 2.0),
                },
                this.orientation,
                { x: 0, y: 0, z: 0 }
            );

            normals.push(newNormal);
        }
        return {
            vertexArr: vertices,
            normalArr: normals,
            faceArr: indices,
            color: color,
        };
    }
}

function triangulateArrow(
    vector: [number, number, number],
    base_radius: number,
    head_radius: number,
    head_length: number,
    resolution: number = 20
): { vertices: XYZ[]; indices: number[] } {
    const [x, y, z] = vector;
    const tip: XYZ = { x, y, z };
    const v_len = Math.sqrt(x * x + y * y + z * z);
    const base_tip: XYZ = {
        x: tip.x * (1 - head_length / v_len),
        y: tip.y * (1 - head_length / v_len),
        z: tip.z * (1 - head_length / v_len),
    };

    // generates a unit circle oriented in the right direction
    const n_vec: XYZ = multXYZ(tip, 1.0 / v_len);

    // Generate an arbitrary vector not collinear with n
    let vx: XYZ;
    if (n_vec.x !== 0.0 || n_vec.y !== 0.0) {
        vx = { x: 0, y: 0, z: 1 };
    } else {
        vx = { x: 0, y: 1, z: 0 };
    }

    // generate orthogonal vectors in the plane defined by nvec
    let u: XYZ = addXYZ(vx, multXYZ(n_vec, -dotXYZ(vx, n_vec)));
    u = multXYZ(u, 1.0 / Math.sqrt(dotXYZ(u, u)));
    const v: XYZ = crossXYZ(u, n_vec);

    // generate n_points in the plane defined by nvec, centered at vec
    const circle_points: XYZ[] = [];
    for (let i = 0; i < resolution; i++) {
        circle_points.push(
            addXYZ(
                multXYZ(u, Math.cos((i * 2 * Math.PI) / resolution)),
                multXYZ(v, Math.sin((i * 2 * Math.PI) / resolution))
            )
        );
    }

    let indices: number[] = [];
    const vertices: XYZ[] = [];

    vertices.push({ x: 0, y: 0, z: 0 });
    // the arrow is built as a surface of revolution, by stacking _|\ motifs
    for (let i = 0; i < resolution; i++) {
        // nb replicated points are needed to get sharp edges
        vertices.push(multXYZ(circle_points[i], base_radius));
        vertices.push(multXYZ(circle_points[i], base_radius));
        vertices.push(addXYZ(multXYZ(circle_points[i], base_radius), base_tip));
        vertices.push(addXYZ(multXYZ(circle_points[i], base_radius), base_tip));
        vertices.push(addXYZ(multXYZ(circle_points[i], head_radius), base_tip));
        vertices.push(addXYZ(multXYZ(circle_points[i], head_radius), base_tip));
        vertices.push(tip);
        const i_seg = 1 + i * 7;
        const i_next = 1 + ((i + 1) % resolution) * 7;
        indices = [
            ...indices,
            ...[
                0,
                i_seg,
                i_next, // cylinder base
                i_seg + 1,
                i_seg + 2,
                i_next + 1,
                i_next + 1,
                i_seg + 2,
                i_next + 2, // cylinder side
                i_seg + 3,
                i_seg + 4,
                i_next + 3,
                i_seg + 4,
                i_next + 4,
                i_next + 3, // head ring
                i_seg + 5,
                i_next + 5,
                i_seg + 6, // tip
            ],
        ];
    }
    return {
        vertices: vertices,
        indices: indices,
    };
}

export class Arrow extends Shape {
    public vector: [number, number, number];
    public base_radius: number;
    public head_radius: number;
    public head_length: number;

    constructor(data: Partial<ArrowData>) {
        super(data);
        assert(data.vector);
        this.vector = [
            this.scale * data.vector[0],
            this.scale * data.vector[1],
            this.scale * data.vector[2],
        ];
        this.base_radius = this.scale * (data.base_radius || 0.1);
        this.head_radius = this.scale * (data.head_radius || 0.15);
        this.head_length = this.scale * (data.head_length || 0.2);
    }

    public static validateParameters(parameters: Record<string, unknown>): string {
        if (!('vector' in parameters)) {
            return '"vector" is required for "arrow" shapes';
        }

        if (!Array.isArray(parameters.vector) || parameters.vector.length !== 3) {
            return '"vector" must be an array with 3 elements for "vector" shapes';
        }

        const [ax, ay, az] = parameters.vector as unknown[];
        if (typeof ax !== 'number' || typeof ay !== 'number' || typeof az !== 'number') {
            return '"vector" elements must be numbers for "vector" shapes';
        }

        return '';
    }

    public outputTo3Dmol(color: $3Dmol.ColorSpec, resolution: number = 20): $3Dmol.CustomShapeSpec {
        const triangulation = triangulateArrow(
            this.vector,
            this.base_radius,
            this.head_radius,
            this.head_length,
            resolution
        );
        const rawVertices = triangulation.vertices;
        const indices = triangulation.indices;
        const vertices: XYZ[] = [];
        const simplices: [number, number, number][] = [];

        for (const v of rawVertices) {
            const newVertex: XYZ = addXYZ(v, this.position);
            vertices.push(newVertex);
        }

        for (let i = 0; i < indices.length / 3; i++) {
            simplices.push([indices[3 * i], indices[3 * i + 1], indices[3 * i + 2]]);
        }

        return {
            vertexArr: vertices,
            normalArr: determineNormals(vertices, simplices),
            faceArr: indices,
            color: color,
        };
    }
}

export class CustomShape extends Shape {
    public vertices: XYZ[];
    public simplices: [number, number, number][];

    constructor(data: Partial<CustomShapeData>) {
        super(data);
        assert(data.vertices);
        assert(data.simplices);

        this.vertices = [];
        for (const v of data.vertices) {
            this.vertices.push({
                x: v[0] * this.scale,
                y: v[1] * this.scale,
                z: v[2] * this.scale,
            });
        }
        this.simplices = data.simplices;
    }

    public static validateParameters(parameters: Record<string, unknown>): string {
        if (!('vertices' in parameters)) {
            return '"vertices" is required for "custom" shapes';
        }

        if (!Array.isArray(parameters.vertices)) {
            return '"vertices" must be an array for "custom" shapes';
        }

        for (const vertex of parameters.vertices as unknown[]) {
            if (!Array.isArray(vertex) || vertex.length !== 3) {
                return 'each entry in "vertices" must be an array with 3 elements for "custom" shapes';
            }

            const [x, y, z] = vertex as unknown[];
            if (typeof x !== 'number' || typeof y !== 'number' || typeof z !== 'number') {
                return 'each entry in "vertices" must be an array of numbers for "custom" shapes';
            }
        }

        if ('simplices' in parameters) {
            if (!Array.isArray(parameters.simplices)) {
                return '"simplices" must be an array for "custom" shapes';
            }

            for (const simplex of parameters.simplices as unknown[]) {
                if (!Array.isArray(simplex) || simplex.length !== 3) {
                    return 'each entry in "simplices" must be an array with 3 elements for "custom" shapes';
                }

                const [x, y, z] = simplex as unknown[];
                if (!isPositiveInteger(x) || !isPositiveInteger(y) || !isPositiveInteger(z)) {
                    return 'each entry in "simplices" must be an array of integers for "custom" shapes';
                }
            }
        }

        if ('orientation' in parameters) {
            const message = validateOrientation(parameters.orientation, 'custom');
            if (message !== '') {
                return message;
            }
        }

        return '';
    }

    public outputTo3Dmol(color: $3Dmol.ColorSpec): $3Dmol.CustomShapeSpec {
        const indices = [];
        const vertices = [];

        for (const v of this.vertices) {
            vertices.push(rotateAndPlace(v, this.orientation, this.position));
        }

        for (const s of this.simplices) {
            for (const ss of s) {
                indices.push(ss);
            }
        }

        return {
            vertexArr: vertices,
            normalArr: determineNormals(vertices, this.simplices),
            faceArr: indices,
            color: color,
        };
    }
}
