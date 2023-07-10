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

import { default as $3Dmol } from '3dmol';
import { Quaternion } from '3dmol';

// Just an interface to enforce XYZ-type coordinates
export interface XYZ {
    x: number;
    y: number;
    z: number;
}

export interface SphereData {
    kind: 'sphere';
    radius: number;
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

export type ShapeData = SphereData | EllipsoidData | CustomShapeData;

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

function rotate_and_place(vertex: XYZ, quaternion: Quaternion, position: XYZ): XYZ {
    const vquat = new Quaternion(vertex.x, vertex.y, vertex.z, 0);
    vquat.multiplyQuaternions(quaternion.clone(), vquat);
    vquat.multiply(quaternion.clone().inverse());

    return addXYZ({ x: vquat.x, y: vquat.y, z: vquat.z }, position);
}

function find_center(a: XYZ[]): XYZ {
    let center: XYZ = { x: 0, y: 0, z: 0 };
    for (const v of a) {
        center = addXYZ(center, v);
    }
    center.x /= a.length;
    center.y /= a.length;
    center.z /= a.length;
    return center;
}

function determine_normals(vertices: XYZ[], simplices: [number, number, number][]): XYZ[] {
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
    constructor(
        pos: [number, number, number],
        orientation: [number, number, number, number] = [0, 0, 0, 1]
    ) {
        const quat_norm = Math.pow(
            Math.pow(orientation[0], 2) +
                Math.pow(orientation[1], 2) +
                Math.pow(orientation[2], 2) +
                Math.pow(orientation[3], 2),
            0.5
        );

        if (Math.abs(quat_norm - 1) > 10e-6) {
            throw Error('Non-normalized quaternions may cause some weird visuals.');
        }

        this.quaternion = new Quaternion(
            orientation[0],
            orientation[1],
            orientation[2],
            orientation[3]
        );
        this.position = { x: pos[0], y: pos[1], z: pos[2] };
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

function triangulate(
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

    //Now triangulation. I believe that it needs a counterclockwise specification of vertices (else I have back-face)
    //each triangle is specified with 3 indices (gl.TRIANGLES).

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

export class Sphere extends Shape {
    public radius: number;

    constructor(pos: [number, number, number] = [0, 0, 0], data: SphereData) {
        super(pos);
        this.radius = data.radius;
    }

    public static validateParams(parameters: Record<string, unknown>): string {
        if ('radius' in parameters) {
            const radius = parameters['radius'];
            if (typeof radius === 'number') {
                return '';
            } else {
                return '`radius` should be a number.';
            }
        } else {
            return '`radius` is a required parameter for type: `sphere`.';
        }
    }
    /**
     ** Following https://computergraphics.stackexchange.com/questions/7426/triangulation-of-vertices-of-an-ellipsoid
     **/
    public outputTo3Dmol(color: $3Dmol.ColorSpec, resolution: number = 20): $3Dmol.CustomShapeSpec {
        const triangulation = triangulate([this.radius, this.radius, this.radius], resolution);
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

    constructor(pos: [number, number, number] = [0, 0, 0], data: EllipsoidData) {
        super(pos, data.orientation);
        this.semiaxes = data.semiaxes;
    }

    public static validateParams(parameters: Record<string, unknown>): string {
        if ('semiaxes' in parameters) {
            const semiaxes = parameters['semiaxes'];
            if (Array.isArray(semiaxes) && semiaxes.length === 3) {
                if (
                    typeof semiaxes[0] === 'number' &&
                    typeof semiaxes[1] === 'number' &&
                    typeof semiaxes[2] === 'number'
                ) {
                    return '';
                } else {
                    return '`semiaxes` must all be numbers.';
                }
            } else {
                return '`semiaxes` must be a length 3 array.';
            }
        } else {
            return '`semiaxes` is a required parameter for type `ellipsoid`.';
        }
    }

    /**
     **
     **/
    public outputTo3Dmol(color: $3Dmol.ColorSpec, resolution: number = 20): $3Dmol.CustomShapeSpec {
        const triangulation = triangulate(this.semiaxes, resolution);
        const rawVertices = triangulation.vertices;
        const indices = triangulation.indices;
        const vertices: XYZ[] = [];
        const normals: XYZ[] = [];

        for (const v of rawVertices) {
            const newVertex: XYZ = rotate_and_place(v, this.quaternion, this.position);
            vertices.push(newVertex);

            const relativeVertex = subXYZ(newVertex, this.position);
            const newNormal: XYZ = {
                x: relativeVertex.x / Math.pow(this.semiaxes[0], 2.0),
                y: relativeVertex.y / Math.pow(this.semiaxes[1], 2.0),
                z: relativeVertex.z / Math.pow(this.semiaxes[2], 2.0),
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

    constructor(pos: [number, number, number] = [0, 0, 0], data: CustomShapeData) {
        super(pos, data.orientation);

        this.vertices = [];
        for (const v of data.vertices) {
            this.vertices.push({ x: v[0], y: v[1], z: v[2] });
        }
        this.simplices = data.simplices;
    }

    public static validateParams(parameters: Record<string, unknown>): string {
        if ('vertices' in parameters) {
            const vertices = parameters['vertices'];
            if (Array.isArray(vertices)) {
                for (const v of vertices) {
                    if (Array.isArray(v) && v.length === 3) {
                        if (
                            typeof v[0] === 'number' &&
                            typeof v[1] === 'number' &&
                            typeof v[2] === 'number'
                        ) {
                            // pass
                        } else {
                            return 'All `vertices` must be numbers.';
                        }
                    } else {
                        return 'Each element of `vertices` a length 3 array.';
                    }
                }
            } else {
                return '`vertices` an array.';
            }
        } else {
            return '`vertices` is a required parameter for type `Poly3D`.';
        }

        return '';
    }

    /**
     **
     **/
    public outputTo3Dmol(color: $3Dmol.ColorSpec): $3Dmol.CustomShapeSpec {
        const indices = [];
        const vertices = [];

        const center = find_center(this.vertices);

        for (const v of this.vertices) {
            vertices.push(rotate_and_place(subXYZ(v, center), this.quaternion, this.position));
        }

        for (const s of this.simplices) {
            for (const ss of s) {
                indices.push(ss);
            }
        }

        return {
            vertexArr: vertices,
            normalArr: determine_normals(vertices, this.simplices),
            faceArr: indices,
            color: color,
        };
    }
}
