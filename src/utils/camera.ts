import assert from 'assert';

/**
 * This module contains utility functions to handle the camera settings
 * for Plotly.js and 3dmol.js viewers, including conversions between their
 * respective camera representations.
 *
 * @packageDocumentation
 * @module utils
 */

import { Vector3D, norm } from '../structure/linalg';
import { Bounds } from '../utils';

export interface CameraState {
    eye: { x: number; y: number; z: number };
    center: { x: number; y: number; z: number };
    up: { x: number; y: number; z: number };
    zoom: number;
}

export interface PlotlyState {
    camera: {
        eye: { x: number; y: number; z: number };
        center: { x: number; y: number; z: number };
        up: { x: number; y: number; z: number };
        projection: { type: string } | undefined;
    };
    aspectratio: { x: number; y: number; z: number };
}

// 3Dmol view: [cx, cy, cz, distance, qx, qy, qz, qw]
export type ViewState = [number, number, number, number, number, number, number, number];

function applyQuat(v: Vector3D, q: [number, number, number, number]): Vector3D {
    const [x, y, z] = v;
    const [qx, qy, qz, qw] = q;

    // calculate quat * vec
    const ix = qw * x + qy * z - qz * y;
    const iy = qw * y + qz * x - qx * z;
    const iz = qw * z + qx * y - qy * x;
    const iw = -qx * x - qy * y - qz * z;

    // calculate result * inverse quat
    return [
        ix * qw + iw * -qx + iy * -qz - iz * -qy,
        iy * qw + iw * -qy + iz * -qx - ix * -qz,
        iz * qw + iw * -qz + ix * -qy - iy * -qx,
    ];
}

/**
 * Converts 3Dmol view array to Plotly-style camera object.
 */
export function viewToCamera(view: ViewState): CameraState {
    const center = { x: view[0], y: view[1], z: view[2] };
    // Observed layout: [cx, cy, cz, distance, x, y, z, w]
    const dist = view[3];
    const q: [number, number, number, number] = [view[4], view[5], view[6], view[7]];

    // Normalize quaternion to avoid scaling issues
    const qLen = Math.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    if (qLen > 0) {
        q[0] /= qLen;
        q[1] /= qLen;
        q[2] /= qLen;
        q[3] /= qLen;
    } else {
        q[3] = 1; // Identity
    }

    // Rotating (0, 1, 0) by q gives up vector.
    const upVec = applyQuat([0, 1, 0], q);

    // Rotating (0, 0, 1) by q gives vector pointing towards eye from center.
    const eyeVec = applyQuat([0, 0, 1], q);

    return {
        center: { ...center },
        eye: {
            x: center.x + eyeVec[0] * dist,
            y: center.y + eyeVec[1] * dist,
            z: center.z + eyeVec[2] * dist,
        },
        up: { x: upVec[0], y: upVec[1], z: upVec[2] },
        zoom: dist,
    };
}

/**
 * Converts Plotly-style camera object to 3Dmol view array.
 */
export function cameraToView(camera: CameraState): ViewState {
    const center = camera.center;
    const eye = camera.eye;
    const up = camera.up;

    // Vector from center to eye
    const zAxis: Vector3D = [eye.x - center.x, eye.y - center.y, eye.z - center.z];
    const zLen = norm(zAxis);
    // Normalize zAxis
    if (zLen > 0) {
        zAxis[0] /= zLen;
        zAxis[1] /= zLen;
        zAxis[2] /= zLen;
    } else {
        zAxis[2] = 1;
    }

    // Up vector (yAxis)
    let yAxis: Vector3D = [up.x, up.y, up.z];
    const yLen = norm(yAxis);
    if (yLen > 0) {
        yAxis[0] /= yLen;
        yAxis[1] /= yLen;
        yAxis[2] /= yLen;
    } else {
        yAxis[1] = 1;
    }

    // X axis = Y cross Z
    const xAxis: Vector3D = [
        yAxis[1] * zAxis[2] - yAxis[2] * zAxis[1],
        yAxis[2] * zAxis[0] - yAxis[0] * zAxis[2],
        yAxis[0] * zAxis[1] - yAxis[1] * zAxis[0],
    ];
    // Recompute Y = Z cross X to ensure orthogonality
    yAxis = [
        zAxis[1] * xAxis[2] - zAxis[2] * xAxis[1],
        zAxis[2] * xAxis[0] - zAxis[0] * xAxis[2],
        zAxis[0] * xAxis[1] - zAxis[1] * xAxis[0],
    ];

    // Rotation Matrix [X Y Z]
    const m00 = xAxis[0],
        m01 = yAxis[0],
        m02 = zAxis[0];
    const m10 = xAxis[1],
        m11 = yAxis[1],
        m12 = zAxis[1];
    const m20 = xAxis[2],
        m21 = yAxis[2],
        m22 = zAxis[2];

    // Convert Matrix to Quaternion
    const trace = m00 + m11 + m22;
    let qw, qx, qy, qz;

    if (trace > 0) {
        const s = 0.5 / Math.sqrt(trace + 1.0);
        qw = 0.25 / s;
        qx = (m21 - m12) * s;
        qy = (m02 - m20) * s;
        qz = (m10 - m01) * s;
    } else {
        if (m00 > m11 && m00 > m22) {
            const s = 2.0 * Math.sqrt(1.0 + m00 - m11 - m22);
            qw = (m21 - m12) / s;
            qx = 0.25 * s;
            qy = (m01 + m10) / s;
            qz = (m02 + m20) / s;
        } else if (m11 > m22) {
            const s = 2.0 * Math.sqrt(1.0 + m11 - m00 - m22);
            qw = (m02 - m20) / s;
            qx = (m01 + m10) / s;
            qy = 0.25 * s;
            qz = (m12 + m21) / s;
        } else {
            const s = 2.0 * Math.sqrt(1.0 + m22 - m00 - m11);
            qw = (m10 - m01) / s;
            qx = (m02 + m20) / s;
            qy = (m12 + m21) / s;
            qz = 0.25 * s;
        }
    }

    return [center.x, center.y, center.z, camera.zoom, qx, qy, qz, qw];
}

/** 4x4 Matrix for 3D projection, row-major */
export type Matrix4 = [
    number,
    number,
    number,
    number,
    number,
    number,
    number,
    number,
    number,
    number,
    number,
    number,
    number,
    number,
    number,
    number,
];

/**
 * Computes the View Matrix from Camera State.
 * (Column-major format for WebGL compatibility, though we calculate row-major logic above)
 * We'll use standard math: M = Translation * Rotation
 */
export function getLookAtMatrix(camera: CameraState): Matrix4 {
    const { eye, center, up } = camera;

    // Z = normalize(eye - center)
    let z0 = center.x - eye.x;
    let z1 = center.y - eye.y;
    let z2 = center.z - eye.z;
    let len = Math.sqrt(z0 * z0 + z1 * z1 + z2 * z2);
    if (len > 0) {
        z0 /= len;
        z1 /= len;
        z2 /= len;
    }

    // X = normalize(cross(up, Z))
    let x0 = up.y * z2 - up.z * z1;
    let x1 = up.z * z0 - up.x * z2;
    let x2 = up.x * z1 - up.y * z0;
    len = Math.sqrt(x0 * x0 + x1 * x1 + x2 * x2);
    if (len > 0) {
        x0 /= len;
        x1 /= len;
        x2 /= len;
    }

    // Y = cross(Z, X)
    const y0 = z1 * x2 - z2 * x1;
    const y1 = z2 * x0 - z0 * x2;
    const y2 = z0 * x1 - z1 * x0;

    // Translation (dot products)
    const dotX = -(x0 * eye.x + x1 * eye.y + x2 * eye.z);
    const dotY = -(y0 * eye.x + y1 * eye.y + y2 * eye.z);
    const dotZ = -(z0 * eye.x + z1 * eye.y + z2 * eye.z);

    // Matrix 4x4
    return [x0, y0, z0, 0, x1, y1, z1, 0, x2, y2, z2, 0, dotX, dotY, dotZ, 1];
}

/**
 * Projects a 3D point to 2D screen space coordinates [-1, 1].
 * Uses orthographic projection based on camera state.
 */
export function projectPoints(
    xValues: number[],
    yValues: number[],
    zValues: number[],
    camera: CameraState,
    bounds: Bounds
): { x: number[]; y: number[] } {
    assert(bounds.z !== undefined);

    const viewMatrix = getLookAtMatrix(camera);
    const zoomFactor = camera.zoom;

    // Get constants to determine the scaled coordinates of the points
    const mx = (bounds.x[0] + bounds.x[1]) / 2;
    const my = (bounds.y[0] + bounds.y[1]) / 2;
    const mz = (bounds.z[0] + bounds.z[1]) / 2;
    const dx = (bounds.x[1] - bounds.x[0]) / 2 / zoomFactor;
    const dy = (bounds.y[1] - bounds.y[0]) / 2 / zoomFactor;
    const dz = (bounds.z[1] - bounds.z[0]) / 2 / zoomFactor;

    const xProj: number[] = [];
    const yProj: number[] = [];

    for (let i = 0; i < xValues.length; i++) {
        // Scale and refer to camera center
        const x = (xValues[i] - mx) / dx - camera.center.x;
        const y = (yValues[i] - my) / dy - camera.center.y;
        const z = (zValues[i] - mz) / dz - camera.center.z;

        // Apply view matrix (vZ is not needed)
        const vX = viewMatrix[0] * x + viewMatrix[4] * y + viewMatrix[8] * z + viewMatrix[12];
        const vY = viewMatrix[1] * x + viewMatrix[5] * y + viewMatrix[9] * z + viewMatrix[13];

        xProj.push(vX);
        yProj.push(vY);
    }

    return { x: xProj, y: yProj };
}

/** Convert internal camera settings to Plotly format */
export function cameraToPlotly(camera: CameraState): PlotlyState {
    const plotlyUpdate: PlotlyState = {
        camera: {
            eye: camera.eye,
            center: camera.center,
            up: camera.up,
            projection: { type: 'orthographic' },
        },
        aspectratio: { x: camera.zoom, y: camera.zoom, z: camera.zoom },
    };

    return plotlyUpdate;
}

/** Convert Plotly scene format to internal settings */
export function plotlyToCamera(plotlyUpdate: PlotlyState): CameraState {
    const camera: CameraState = {
        eye: plotlyUpdate.camera.eye,
        center: plotlyUpdate.camera.center,
        up: plotlyUpdate.camera.up,
        zoom:
            (plotlyUpdate.aspectratio.x + plotlyUpdate.aspectratio.y + plotlyUpdate.aspectratio.z) /
            3,
    };

    return camera;
}

/**
 * Validates that the input object matches the CameraState interface.
 * Throws an error if the validation fails.
 */
export function validateCamera(camera: CameraState): void {
    const data = camera as unknown as Record<string, unknown>;

    if (typeof data !== 'object' || data === null) {
        throw Error('invalid type for camera, expected object');
    }

    // Check top level keys
    for (const key of ['eye', 'center', 'up']) {
        if (!(key in data)) {
            throw Error(`missing key '${key}' in camera`);
        }
        const vec = data[key] as Record<string, unknown>;
        if (typeof vec !== 'object' || vec === null) {
            throw Error(`invalid type for camera.${key}, expected object`);
        }
        for (const subkey of ['x', 'y', 'z']) {
            if (!(subkey in vec)) {
                throw Error(`missing key '${subkey}' in camera.${key}`);
            }
            if (typeof vec[subkey] !== 'number') {
                throw Error(`invalid type for camera.${key}.${subkey}, expected number`);
            }
        }
    }

    if (!('zoom' in data)) {
        throw Error("missing key 'zoom' in camera");
    }
    if (typeof data.zoom !== 'number') {
        throw Error('invalid type for camera.zoom, expected number');
    }
}
