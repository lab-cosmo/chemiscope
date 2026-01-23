/**
 * @packageDocumentation
 * @module map
 */

import { Bounds, arrayMaxMin } from '../utils';
import { CameraState, projectPoints } from '../utils/camera';

function computeBounds(xValues: number[], yValues: number[], zValues: number[] | null): Bounds {
    const xRange = arrayMaxMin(xValues);
    const yRange = arrayMaxMin(yValues);
    const zRange = zValues ? arrayMaxMin(zValues) : { min: 0, max: 1 };

    return {
        x: [xRange.min, xRange.max],
        y: [yRange.min, yRange.max],
        z: [zRange.min, zRange.max],
    };
}

function extractGridIds(grid: Int32Array): number[] {
    const indices: number[] = [];
    for (let i = 0; i < grid.length; i++) {
        if (grid[i] !== -1) {
            indices.push(grid[i]);
        }
    }
    return indices.sort((a, b) => a - b);
}

function binPointsByDataCoordinates(
    xValues: number[],
    yValues: number[],
    zValues: number[] | null,
    maxPoints: number,
    bounds?: Bounds
): number[] {
    const is3D = zValues !== null;
    const clipToBounds = bounds !== undefined;
    const dataBounds = clipToBounds ? bounds : computeBounds(xValues, yValues, zValues);

    // Filter points within bounds
    const visibleIds: number[] = [];
    for (let i = 0; i < xValues.length; i++) {
        const x = xValues[i];
        const y = yValues[i];

        if (clipToBounds) {
            const [xMin, xMax] = dataBounds.x;
            const [yMin, yMax] = dataBounds.y;

            if (x < xMin || x > xMax || y < yMin || y > yMax) {
                continue;
            }

            if (is3D && zValues) {
                const z = zValues[i];
                const [zMin, zMax] = dataBounds.z ?? [0, 1];
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
    const [zMin, zMax] = dataBounds.z ?? [0, 1];

    const xRange = xMax - xMin || 1;
    const yRange = yMax - yMin || 1;
    const zRange = zMax - zMin || 1;

    const gridResolution = is3D ? Math.ceil(Math.cbrt(maxPoints)) : Math.ceil(Math.sqrt(maxPoints));
    const gridSize = is3D ? gridResolution ** 3 : gridResolution ** 2;
    const grid = new Int32Array(gridSize).fill(-1);

    // Pre-calculate factors
    const xScale = gridResolution / xRange;
    const yScale = gridResolution / yRange;
    const zScale = zRange === 0 ? 0 : gridResolution / zRange;

    const clamp = (i: number, size: number) => Math.max(0, Math.min(i, size - 1));

    const strideY = gridResolution;
    const strideZ = gridResolution * gridResolution;

    for (const id of visibleIds) {
        const xi = clamp(Math.floor((xValues[id] - xMin) * xScale), gridResolution);
        const yi = clamp(Math.floor((yValues[id] - yMin) * yScale), gridResolution);

        let cellId = xi + yi * strideY;
        if (is3D && zValues) {
            const zi = clamp(Math.floor((zValues[id] - zMin) * zScale), gridResolution);
            cellId += zi * strideZ;
        }

        if (grid[cellId] === -1) {
            grid[cellId] = id;
        }
    }

    return extractGridIds(grid);
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
    const gridSize = clipSize * 2;
    const invStep = gridResolution / gridSize;

    for (const id of visibleIds) {
        const x = projections.x[id];
        const y = projections.y[id];

        // Convert screen coordinates to grid indices
        const xi = Math.floor((x + clipSize) * invStep);
        const yi = Math.floor((y + clipSize) * invStep);

        if (xi < 0 || xi >= gridResolution || yi < 0 || yi >= gridResolution) {
            continue;
        }

        const cellId = xi + yi * gridResolution;
        if (grid[cellId] === -1) {
            grid[cellId] = id;
        }
    }

    return extractGridIds(grid);
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

    // Coarse pass
    const coarseSpatialIds = binPointsByDataCoordinates(
        xValues,
        yValues,
        zValues,
        Math.floor(maxPoints / 10)
    );

    // Fine pass
    let fineDetailIds: number[];
    if (is3D && zValues && camera && bounds) {
        const effectiveBounds = bounds ?? computeBounds(xValues, yValues, zValues);
        fineDetailIds = binPointsByCameraProjection(
            xValues,
            yValues,
            zValues,
            camera,
            effectiveBounds,
            Math.floor(maxPoints / 2)
        );
    } else {
        fineDetailIds = binPointsByDataCoordinates(xValues, yValues, zValues, maxPoints, bounds);
    }

    // Combine and delete duplicates
    return [...new Set([...coarseSpatialIds, ...fineDetailIds])].sort((a, b) => a - b);
}
