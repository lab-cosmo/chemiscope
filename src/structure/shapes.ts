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

// Just an interface to enforce XYZ-type coordinates
export interface XYZ {
    x: number;
    y: number;
    z: number;
}

export interface BaseShapeData {
    position: [number, number, number];    
    orientation: [number, number, number, number];
    color: $3Dmol.ColorSpec;
}

export interface BaseShapeParameters<T> {
    kind: string;
    parameters: {
        global: Partial<T>,
        structure?: Array<Partial<T>>,
        atom?: Array<Partial<T>>,
    }
}

export interface SphereData extends BaseShapeData{
    radius: number;
}

export interface SphereParameters  extends BaseShapeParameters<SphereData> {
    kind: 'sphere';
}

// Interface for ellipsoidal data, where
// orientation is stored in the (w, x, y, z) convention
export interface EllipsoidData {
    kind: 'ellipsoid';
    semiaxes: [number, number, number];
    orientation?: [number, number, number, number];
}

// Interface for polytope data, where
// orientation is stored in the (w, x, y, z) convention
// and simplices refers to the indices of the facets
export interface CustomShapeData {
    kind: 'custom';
    vertices: [number, number, number][];
    simplices: [number, number, number][];
    orientation?: [number, number, number, number];
}

export type ShapeData = SphereData ; //| EllipsoidData | CustomShapeData;
export type ShapeParameters = SphereParameters ; 

function addXYZ(a: XYZ, b: XYZ): XYZ {
    return { x: a.x + b.x, y: a.y + b.y, z: a.z + b.z };
}

function subXYZ(a: XYZ, b: XYZ): XYZ {
    return { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z };
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
    public quaternion: Quaternion;

    // orientation is passed to 3dmol in the (x, y, z, w) convention
    constructor( data: Partial<ShapeData> ) {
        let [qx, qy, qz, qw] = data.orientation || [0,0,0,1];
        this.quaternion = new Quaternion(qx, qy, qz, qw);
        let [x, y, z] = data.position || [0,0,0];
        this.position = { x, y, z };
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
        super({"position": data.position});
        this.radius = data.radius || 1.0 ;        
    }

    public static validateParameters(parameters: Record<string, unknown>): string {
        assert(parameters.kind === 'sphere');

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

    constructor(position: [number, number, number] = [0, 0, 0], data: EllipsoidData) {
        super({});
        //super(position, data.orientation);
        this.semiaxes = data.semiaxes;
    }

    public static validateParameters(parameters: Record<string, unknown>): string {
        assert(parameters.kind === 'ellipsoid');

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
            const newVertex: XYZ = rotateAndPlace(v, this.quaternion, this.position);
            vertices.push(newVertex);

            const relativeVertex = subXYZ(newVertex, this.position);
            const newNormal: XYZ = {
                /*
                I think these are wrong because the vertex has now been rotated so it's not aligned with the axes
                x: relativeVertex.x / Math.pow(this.semiaxes[0], 2.0),
                y: relativeVertex.y / Math.pow(this.semiaxes[1], 2.0),
                z: relativeVertex.z / Math.pow(this.semiaxes[2], 2.0),
                */
                // this is also wrong, but possibly less wrong. 
                x: relativeVertex.x,
                y: relativeVertex.y,
                z: relativeVertex.z,
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

export class CustomShape extends Shape {
    public vertices: XYZ[];
    public simplices: [number, number, number][];

    constructor(position: [number, number, number] = [0, 0, 0], data: CustomShapeData) {
        super({}); //position, data.orientation);

        this.vertices = [];
        for (const v of data.vertices) {
            this.vertices.push({ x: v[0], y: v[1], z: v[2] });
        }
        this.simplices = data.simplices;
    }

    public static validateParameters(parameters: Record<string, unknown>): string {
        assert(parameters.kind === 'custom');

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
            vertices.push(rotateAndPlace(v, this.quaternion, this.position));
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
