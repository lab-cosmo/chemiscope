/**
 * @packageDocumentation
 * @module map
 */

import { Bounds, arrayMaxMin } from '../utils';
import { CameraState, projectPoints } from '../utils/camera';

/** Clamps a value to [0, max-1] */
function clampToGrid(value: number, max: number): number {
    return Math.max(0, Math.min(max - 1, value));
}

/** Computes bounds from data if not provided */
function getOrComputeBounds(
    xValues: number[],
    yValues: number[],
    zValues: number[] | null,
    bounds?: Bounds
): Bounds {
    if (bounds) {
        return bounds;
    }

    const xRange = arrayMaxMin(xValues);
    const yRange = arrayMaxMin(yValues);
    const zRange = zValues ? arrayMaxMin(zValues) : { min: 0, max: 1 };

    return {
        x: [xRange.min, xRange.max],
        y: [yRange.min, yRange.max],
        z: [zRange.min, zRange.max],
    };
}

/** Extracts filled grid cells into a sorted index array */
function extractGridIndices(grid: Int32Array): number[] {
    const indices: number[] = [];
    for (let i = 0; i < grid.length; i++) {
        if (grid[i] !== -1) {
            indices.push(grid[i]);
        }
    }
    return indices.sort((a, b) => a - b);
}

/** Bins points based on their data coordinates (X, Y, Z) */
function binPointsByDataCoordinates(
    xValues: number[],
    yValues: number[],
    zValues: number[] | null,
    maxPoints: number,
    bounds?: Bounds
): number[] {
    const is3D = zValues !== null;
    const dataBounds = getOrComputeBounds(xValues, yValues, zValues, bounds);

    // Filter points within bounds
    const visibleIds: number[] = [];
    for (let i = 0; i < xValues.length; i++) {
        const x = xValues[i];
        const y = yValues[i];

        if (bounds) {
            const [xMin, xMax] = dataBounds.x;
            const [yMin, yMax] = dataBounds.y;

            if (x < xMin || x > xMax || y < yMin || y > yMax) {
                continue;
            }

            if (is3D && zValues && dataBounds.z) {
                const z = zValues[i];
                const [zMin, zMax] = dataBounds.z;
                if (z < zMin || z > zMax) {
                    continue;
                }
            }
        }
        visibleIds.push(i);
    }

    if (visibleIds.length < maxPoints) {
        return visibleIds;
    }

    // Create spatial grid
    const [xMin, xMax] = dataBounds.x;
    const [yMin, yMax] = dataBounds.y;
    const [zMin, zMax] = dataBounds.z || [0, 1];

    const xRange = xMax - xMin || 1;
    const yRange = yMax - yMin || 1;
    const zRange = zMax - zMin || 1;

    const gridResolution = is3D ? Math.ceil(Math.cbrt(maxPoints)) : Math.ceil(Math.sqrt(maxPoints));
    const gridSize = is3D ? gridResolution ** 3 : gridResolution ** 2;
    const grid = new Int32Array(gridSize).fill(-1);

    // Pre-calculate factors
    const xScale = gridResolution / xRange;
    const yScale = gridResolution / yRange;
    const zScale = gridResolution / zRange;

    for (const id of visibleIds) {
        const x = xValues[id];
        const y = yValues[id];

        const xi = clampToGrid(Math.floor((x - xMin) * xScale), gridResolution);
        const yi = clampToGrid(Math.floor((y - yMin) * yScale), gridResolution);

        let cellIndex: number;
        if (is3D && zValues) {
            const z = zValues[id];
            const zi = clampToGrid(Math.floor((z - zMin) * zScale), gridResolution);
            cellIndex = xi + yi * gridResolution + zi * gridResolution ** 2;
        } else {
            cellIndex = xi + yi * gridResolution;
        }

        if (grid[cellIndex] === -1) {
            grid[cellIndex] = id;
        }
    }

    return extractGridIndices(grid);
}

/** Bins points based on camera projection (3D only) */
function binPointsByCameraProjection(
    xValues: number[],
    yValues: number[],
    zValues: number[],
    camera: CameraState,
    bounds: Bounds,
    maxPoints: number
): number[] {
    // Project all 3D points to 2D coordinates
    const projections = projectPoints(xValues, yValues, zValues, camera, bounds);

    const viewClipPadding = 1.1;
    const clipSize = 2.0 * viewClipPadding;

    // Find points within the camera view
    const visibleIds: number[] = [];
    for (let i = 0; i < xValues.length; i++) {
        const x = projections.x[i];
        const y = projections.y[i];

        if (Math.abs(x) <= clipSize && Math.abs(y) <= clipSize) {
            visibleIds.push(i);
        }
    }

    if (visibleIds.length < maxPoints) {
        return visibleIds;
    }

    // Create 2D grid in screen space
    const gridResolution = Math.ceil(Math.sqrt(maxPoints));
    const grid = new Int32Array(gridResolution * gridResolution).fill(-1);
    const gridStep = (clipSize * 2) / gridResolution;

    for (const id of visibleIds) {
        const x = projections.x[id];
        const y = projections.y[id];

        // Convert screen coordinates to grid indices
        const xi = Math.floor((x + clipSize) / gridStep);
        const yi = Math.floor((y + clipSize) / gridStep);

        if (xi < 0 || xi >= gridResolution || yi < 0 || yi >= gridResolution) {
            continue;
        }

        const cellIndex = xi + yi * gridResolution;
        if (grid[cellIndex] === -1) {
            grid[cellIndex] = id;
        }
    }

    return extractGridIndices(grid);
}

/**
 * Computes the subset of points to display with LOD
 *
 * @param xValues Array of X coordinates
 * @param yValues Array of Y coordinates
 * @param zValues Array of Z coordinates (or null for 2D)
 * @param maxPoints Maximum number of points to display
 * @param bounds Optional boundaries to clip the data
 * @param camera Optional camera state for 3D screen-space LOD
 * @returns Sorted array of indices to display
 */
export function computeLODIndices(
    xValues: number[],
    yValues: number[],
    zValues: number[] | null,
    maxPoints: number,
    bounds?: Bounds,
    camera?: CameraState
): number[] | null {
    if (xValues.length <= maxPoints) {
        return null;
    }

    const is3D = zValues !== null;

    const coarseSpatialIds = binPointsByDataCoordinates(
        xValues,
        yValues,
        zValues,
        Math.floor(maxPoints / 10)
    );

    let fineDetailIds: number[];
    if (is3D && zValues && camera && bounds) {
        // Use camera projection binning
        const effectiveBounds = getOrComputeBounds(xValues, yValues, zValues, bounds);
        fineDetailIds = binPointsByCameraProjection(
            xValues,
            yValues,
            zValues,
            camera,
            effectiveBounds,
            Math.floor(maxPoints / 2)
        );
    } else {
        // Use coordinate from data with current zoom bounds
        fineDetailIds = binPointsByDataCoordinates(xValues, yValues, zValues, maxPoints, bounds);
    }

    // Combine and delete duplicates
    return [...new Set([...coarseSpatialIds, ...fineDetailIds])].sort((a, b) => a - b);
}
