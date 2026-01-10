/**
 * This module contains utility functions to handle the camera settings
 * for Plotly.js and 3dmol.js viewers, including conversions between their
 * respective camera representations.
 * 
 * @packageDocumentation
 * @module utils
 */

import { Vector3D, norm } from '../structure/linalg';

export interface CameraState {
    eye: { x: number; y: number; z: number };
    center: { x: number; y: number; z: number };
    up: { x: number; y: number; z: number };
    zoom: number;
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
    const m00 = xAxis[0], m01 = yAxis[0], m02 = zAxis[0];
    const m10 = xAxis[1], m11 = yAxis[1], m12 = zAxis[1];
    const m20 = xAxis[2], m21 = yAxis[2], m22 = zAxis[2];

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

    return [center.x, center.y, center.z, zLen, qx, qy, qz, qw];
}

/** Convert internal camera settings to Plotly format */
// eslint-disable-next-line @typescript-eslint/no-explicit-any
export function cameraToPlotly(camera: any): any {
    const plotlyCamera = { ...camera };
    if (plotlyCamera.zoom !== undefined) {
        plotlyCamera.projection = { type: 'orthographic', scale: plotlyCamera.zoom };
        delete plotlyCamera.zoom;
    } else {
        plotlyCamera.projection = { type: 'orthographic' };
    }
    return plotlyCamera;
}

/** Convert Plotly camera format to internal settings */
// eslint-disable-next-line @typescript-eslint/no-explicit-any
export function plotlyToCamera(plotlyCamera: any, scale?: number): any {
    const camera = { ...plotlyCamera };

    // Use provided scale or try to find it in projection
    let zoom = scale;
    if (zoom === undefined && camera.projection && camera.projection.scale !== undefined) {
        zoom = camera.projection.scale;
    }

    if (zoom !== undefined) {
        camera.zoom = zoom;
    }

    delete camera.projection;
    return camera;
}
