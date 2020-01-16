/**
 * Some very basic linear algebra functionalities for 3D vectors and 3x3
 * matrixes.
 *
 * @packageDocumentation
 * @module structure.linalg
 */

export type Vector3D = [number, number, number];
export type Matrix = [Vector3D, Vector3D, Vector3D];

export function determinant(matrix: Matrix): number {
    let determinant = 0.0;
    determinant += matrix[0][0] * matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2];
    determinant -= matrix[0][1] * matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0];
    determinant += matrix[0][2] * matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0];
    return determinant;
}

export function invert(matrix: Matrix): Matrix {
    let invdet = 1.0 / determinant(matrix);
    let xx = (matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2]) * invdet;
    let xy = (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) * invdet;
    let xz = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) * invdet;

    let yx = (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) * invdet;
    let yy = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) * invdet;
    let yz = (matrix[1][0] * matrix[0][2] - matrix[0][0] * matrix[1][2]) * invdet;

    let zx = (matrix[1][0] * matrix[2][1] - matrix[2][0] * matrix[1][1]) * invdet;
    let zy = (matrix[2][0] * matrix[0][1] - matrix[0][0] * matrix[2][1]) * invdet;
    let zz = (matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]) * invdet;

    return [
        [xx, xy, xz],
        [yx, yy, yz],
        [zx, zy, zz],
    ];
}

export function dot(matrix: Matrix, vector: Vector3D): Vector3D {
    return [
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    ];
}

function vec_dot(lhs: Vector3D, rhs: Vector3D): number {
    return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

export function angle(lhs: Vector3D, rhs: Vector3D) {
    const PI = 3.141592653589793;
    return Math.acos(vec_dot(lhs, rhs) / (norm(lhs) * norm(rhs))) * 180. / PI
}

export function norm(vector: Vector3D): number {
    return Math.sqrt(vec_dot(vector, vector));
}
