/** @ignore */
/**
 * Some very basic linear algebra functionalities for 3D vectors and 3x3
 * matrixes.
 */

import assert from 'assert';

export type Vector3D = [number, number, number];
export type Matrix = [Vector3D, Vector3D, Vector3D];

export function determinant(matrix: Matrix): number {
    let result = 0.0;
    result += matrix[0][0] * matrix[1][1] * matrix[2][2];
    result += matrix[0][1] * matrix[1][2] * matrix[2][0];
    result += matrix[0][2] * matrix[1][0] * matrix[2][1];
    result -= matrix[0][0] * matrix[2][1] * matrix[1][2];
    result -= matrix[0][1] * matrix[1][0] * matrix[2][2];
    result -= matrix[0][2] * matrix[1][1] * matrix[2][0];
    return result;
}

export function invert(matrix: Matrix): Matrix {
    const det = determinant(matrix);
    assert(Math.abs(det) > 1e-12);
    const xx = (matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2]) / det;
    const xy = (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / det;
    const xz = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / det;

    const yx = (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / det;
    const yy = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / det;
    const yz = (matrix[1][0] * matrix[0][2] - matrix[0][0] * matrix[1][2]) / det;

    const zx = (matrix[1][0] * matrix[2][1] - matrix[2][0] * matrix[1][1]) / det;
    const zy = (matrix[2][0] * matrix[0][1] - matrix[0][0] * matrix[2][1]) / det;
    const zz = (matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]) / det;

    return [
        [xx, xy, xz],
        [yx, yy, yz],
        [zx, zy, zz],
    ];
}

export function dot(matrix: Matrix, vector: Vector3D): Vector3D {
    return [
        matrix[0][0] * vector[0] + matrix[1][0] * vector[1] + matrix[2][0] * vector[2],
        matrix[0][1] * vector[0] + matrix[1][1] * vector[1] + matrix[2][1] * vector[2],
        matrix[0][2] * vector[0] + matrix[1][2] * vector[1] + matrix[2][2] * vector[2],
    ];
}

function vec_dot(lhs: Vector3D, rhs: Vector3D): number {
    return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

export function angle(lhs: Vector3D, rhs: Vector3D): number {
    const PI = 3.141592653589793;
    return Math.acos(vec_dot(lhs, rhs) / (norm(lhs) * norm(rhs))) * 180 / PI;
}

export function norm(vector: Vector3D): number {
    return Math.sqrt(vec_dot(vector, vector));
}
