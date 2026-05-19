/**
 * @packageDocumentation
 * @module map
 */

import { Bounds, arrayMaxMin } from '../utils';
import { CameraState, projectPoints } from '../utils/camera';

/**
 * Computes LOD based on 2D screen-space projection from 3D
 *
 * @param xValues Array of X coordinates
 * @param yValues Array of Y coordinates
 * @param zValues Array of Z coordinates
 * @param camera Current camera state (eye, center, up, zoom)
 * @param bounds Optional boundaries to clip the data
 * @param maxPoints Maximum number of points to display
 * @param aspectRatio Width/height ratio of plot area
 * @returns Array of indices to display
 */
export function computeScreenSpaceLOD(
    xValues: number[],
    yValues: number[],
    zValues: number[],
    camera: CameraState,
    bounds?: Bounds,
    maxPoints: number = 50000,
    aspectRatio: number = 1.0
): number[] {
    if (bounds === undefined) {
        const xRange = arrayMaxMin(xValues);
        const yRange = arrayMaxMin(yValues);
        const zRange = arrayMaxMin(zValues);
        bounds = {
            x: [xRange.min, xRange.max],
            y: [yRange.min, yRange.max],
            z: [zRange.min, zRange.max],
        };
    }

    const projections = projectPoints(xValues, yValues, zValues, camera, bounds);

    // Bounds slightly larger to avoid popping at edges
    const BASE_CLIP_SIZE = 2.2;
    // Adjust clip sizes based on aspect ratio
    const CLIP_SIZE_X = BASE_CLIP_SIZE * Math.max(1, aspectRatio);
    const CLIP_SIZE_Y = BASE_CLIP_SIZE * Math.max(1, 1 / aspectRatio);

    // Filter visible points
    const visibleIds: number[] = [];
    for (let i = 0; i < xValues.length; i++) {
        if (
            Math.abs(projections.x[i]) <= CLIP_SIZE_X &&
            Math.abs(projections.y[i]) <= CLIP_SIZE_Y
        ) {
            visibleIds.push(i);
        }
    }

    if (visibleIds.length <= maxPoints) {
        return visibleIds;
    }

    // Binning
    const baseBins = Math.ceil(Math.sqrt(maxPoints));
    const binsX = Math.ceil(baseBins * Math.sqrt(aspectRatio));
    const binsY = Math.ceil(baseBins / Math.sqrt(aspectRatio));
    const grid = new Int32Array(binsX * binsY).fill(-1);
    // Per-bin best hash: keep the id with the smallest hash so each cell's
    // pick is a deterministic, evenly-distributed permutation of input order.
    // Selecting the first hit introduces a "first encountered wins" bias
    // that creates artificial structure when the input is sorted.
    const bestHash = new Uint32Array(binsX * binsY).fill(0xffffffff);
    const invStepX = binsX / (CLIP_SIZE_X * 2);
    const invStepY = binsY / (CLIP_SIZE_Y * 2);

    for (const id of visibleIds) {
        const ui = Math.floor((projections.x[id] + CLIP_SIZE_X) * invStepX);
        const vi = Math.floor((projections.y[id] + CLIP_SIZE_Y) * invStepY);

        if (ui < 0 || ui >= binsX || vi < 0 || vi >= binsY) {
            continue;
        }

        const idx = ui + vi * binsX;
        const h = Math.imul(id, 2654435761) >>> 0;
        if (h < bestHash[idx]) {
            bestHash[idx] = h;
            grid[idx] = id;
        }
    }

    const result: number[] = [];
    for (let i = 0; i < grid.length; i++) {
        if (grid[i] !== -1) {
            result.push(grid[i]);
        }
    }
    result.sort((a, b) => a - b);
    return result;
}

/**
 * Computes LOD based on spatial grid binning. Used for 2D plots or to get a global
 * 3D subsampling that does not depend on the camera viewpoint
 *
 * @param xValues Array of X coordinates
 * @param yValues Array of Y coordinates
 * @param zValues Array of Z coordinates (or null for 2D)
 * @param bounds Optional boundaries to clip the data
 * @param maxPoints Maximum number of points to display
 * @returns Array of indices to display
 */
export function computeLODIndices(
    xValues: number[],
    yValues: number[],
    zValues: number[] | null,
    bounds?: Bounds,
    maxPoints: number = 50000
): number[] {
    const is3D = zValues !== null;

    // Determine the range we are binning over
    let xMin: number, xMax: number, yMin: number, yMax: number;
    let zMin = 0;
    let zMax = 1;

    if (bounds) {
        // DYNAMIC: Use the current zoom level provided by bounds
        [xMin, xMax] = bounds.x;
        [yMin, yMax] = bounds.y;
        if (is3D && bounds.z) {
            [zMin, zMax] = bounds.z;
        }
    } else {
        // STATIC: Use the full data range (calculate from data)
        const xRange = arrayMaxMin(xValues);
        xMin = xRange.min;
        xMax = xRange.max;

        const yRange = arrayMaxMin(yValues);
        yMin = yRange.min;
        yMax = yRange.max;

        if (is3D && zValues) {
            const zRange = arrayMaxMin(zValues);
            zMin = zRange.min;
            zMax = zRange.max;
        }
    }

    // Filter visible points
    const visibleIds: number[] = [];
    for (let i = 0; i < xValues.length; i++) {
        const x = xValues[i];
        const y = yValues[i];

        if (bounds) {
            if (x < xMin || x > xMax || y < yMin || y > yMax) {
                continue;
            }

            if (is3D && zValues) {
                const z = zValues[i];
                if (z < zMin || z > zMax) {
                    continue;
                }
            }
        }

        visibleIds.push(i);
    }

    if (visibleIds.length <= maxPoints) {
        return visibleIds;
    }

    // Avoid division by zero
    const xRange = xMax - xMin || 1;
    const yRange = yMax - yMin || 1;
    const zRange = zMax - zMin || 1;

    const bins = is3D ? Math.ceil(Math.cbrt(maxPoints)) : Math.ceil(Math.sqrt(maxPoints));
    const nCells = is3D ? bins ** 3 : bins ** 2;
    const grid = new Int32Array(nCells).fill(-1);
    // Per-bin best hash: keep the id with the smallest hash so each cell's
    // pick is a deterministic, evenly-distributed permutation of input order.
    // This removes the "first encountered wins" bias that creates artificial
    // structure when the input is sorted along a coordinate.
    const bestHash = new Uint32Array(nCells).fill(0xffffffff);

    // Pre-calculate factors
    const xFactor = bins / xRange;
    const yFactor = bins / yRange;
    const zFactor = bins / zRange;

    const clamp = (v: number, max: number) => Math.max(0, Math.min(v, max - 1));

    for (const id of visibleIds) {
        const xi = clamp(Math.floor((xValues[id] - xMin) * xFactor), bins);
        const yi = clamp(Math.floor((yValues[id] - yMin) * yFactor), bins);

        let idx = xi + yi * bins;
        if (is3D && zValues) {
            const zi = clamp(Math.floor((zValues[id] - zMin) * zFactor), bins);
            idx += zi * bins * bins;
        }

        const h = Math.imul(id, 2654435761) >>> 0;
        if (h < bestHash[idx]) {
            bestHash[idx] = h;
            grid[idx] = id;
        }
    }

    const result: number[] = [];
    for (let i = 0; i < grid.length; i++) {
        if (grid[i] !== -1) {
            result.push(grid[i]);
        }
    }
    result.sort((a, b) => a - b);
    return result;
}
