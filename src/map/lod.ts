/**
 * @packageDocumentation
 * @module map
 */

import { arrayMaxMin } from '../utils';
import { CameraState, projectPoints } from '../utils/camera';

/**
 * Computes the subset of points to display based on 2D screen-space projection.
 * This is used for 3D plots where the visible density depends on camera rotation/zoom.
 *
 * @param xValues Array of X coordinates
 * @param yValues Array of Y coordinates
 * @param zValues Array of Z coordinates
 * @param camera  Current camera state (eye, center, up, zoom)
 * @param bounds  Optional boundaries to clip the data (Zoom level)
 * @param threshLOD Maximum number of points to display (default: 50000)
 * @returns Sorted array of indices to display
 */
export function computeScreenSpaceLOD(
    xValues: number[],
    yValues: number[],
    zValues: number[],
    camera: CameraState,
    bounds?: { x: [number, number]; y: [number, number]; z?: [number, number] },
    threshLOD: number = 50000
): number[] {
    if (bounds === undefined) {
        // STATIC: Use the full data range (calculate from data)
        const xRange = arrayMaxMin(xValues);
        const yRange = arrayMaxMin(yValues);
        const zRange = arrayMaxMin(zValues);

        bounds = {
            x: [xRange.min, xRange.max],
            y: [yRange.min, yRange.max],
            z: [zRange.min, zRange.max],
        };
    }

    const nPoints = xValues.length;
    // Calculate View Matrix once
    const projections = projectPoints(xValues, yValues, zValues, camera, bounds);

    // this is the range we use to bin points
    const viewSize = 2.0;
    // Use a slightly larger clip bound to avoid popping at edges.
    // We use half-width (radius) for comparison, so 1.1 covers [-1.1, 1.1]
    const clipSize = viewSize * 1.1;

    // Pass 1: Filter points within the view
    const visibleIndices: number[] = [];
    for (let i = 0; i < nPoints; i++) {
        const u = projections.x[i];
        const v = projections.y[i];
        // Clip points well outside the view
        if (Math.abs(u) <= clipSize && Math.abs(v) <= clipSize) {
            visibleIndices.push(i);
        }
    }

    // Early exit if we have fewer points than the threshold
    if (visibleIndices.length < threshLOD) {
        return visibleIndices;
    }

    // Pass 2: Binning (only for visible points)
    // Grid resolution: we want roughly threshLOD points on screen.
    const bins = Math.ceil(Math.sqrt(threshLOD));
    // Use Int32Array for performance (avoid string allocation)
    // Initialize with -1 to indicate empty cells
    const grid = new Int32Array(bins * bins).fill(-1);

    // The grid covers the range [-clipSize, clipSize]
    const gridSize = clipSize * 2;
    const step = gridSize / bins;

    // Pre-calculate inverse step to replace division with multiplication
    const invStep = 1.0 / step;
    // Offset to shift [-clipSize, clipSize] to [0, gridSize]
    const offset = clipSize;

    for (let j = 0; j < visibleIndices.length; j++) {
        const i = visibleIndices[j];
        const u = projections.x[i];
        const v = projections.y[i];

        // Map to grid coordinates [0, bins-1]
        const ui = Math.floor((u + offset) * invStep);
        const vi = Math.floor((v + offset) * invStep);

        // Safety clamp (floating point errors might put 1.1 slightly over)
        if (ui < 0 || ui >= bins || vi < 0 || vi >= bins) {
            continue;
        }

        const index = ui + vi * bins;
        if (grid[index] === -1) {
            grid[index] = i;
        }
    }

    // Filter -1s and sort
    const indices: number[] = [];
    for (let i = 0; i < grid.length; i++) {
        if (grid[i] !== -1) {
            indices.push(grid[i]);
        }
    }
    indices.sort((a, b) => a - b);
    return indices;
}

/**
 * Computes the subset of points to display based on spatial grid binning.
 *
 * @param xValues Array of X coordinates
 * @param yValues Array of Y coordinates
 * @param zValues Array of Z coordinates (or null for 2D)
 * @param bounds  Optional boundaries to clip the data (Zoom level)
 * @param threshLOD  Maximum number of points to display (default: 50000)
 * @returns Sorted array of indices to display
 */
export function computeLODIndices(
    xValues: number[],
    yValues: number[],
    zValues: number[] | null,
    bounds?: { x: [number, number]; y: [number, number]; z?: [number, number] },
    threshLOD: number = 50000
): number[] {
    const nPoints = xValues.length;
    const is3D = zValues !== null;

    // Determine the range we are binning over
    let xMin: number, xMax: number, yMin: number, yMax: number;
    let zMin = 0;
    let zMax = 0;

    if (bounds) {
        // DYNAMIC: Use the current zoom level provided by bounds
        [xMin, xMax] = bounds.x;
        [yMin, yMax] = bounds.y;
        if (is3D && bounds.z) {
            [zMin, zMax] = bounds.z;
        } else {
            zMin = 0;
            zMax = 1;
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
        } else {
            zMin = 0;
            zMax = 1;
        }
    }

    // Pass 1: Filter points within bounds
    const visibleIndices: number[] = [];
    for (let i = 0; i < nPoints; i++) {
        const xVal = xValues[i];
        const yVal = yValues[i];

        if (bounds) {
            if (xVal < xMin || xVal > xMax || yVal < yMin || yVal > yMax) {
                continue;
            }
            if (is3D && zValues) {
                const zVal = zValues[i];
                if (zVal < zMin || zVal > zMax) {
                    continue;
                }
            }
        }
        visibleIndices.push(i);
    }

    // Early exit
    if (visibleIndices.length < threshLOD) {
        return visibleIndices;
    }

    // Pass 2: Binning (only for clipped points)
    const xRange = xMax - xMin || 1;
    const yRange = yMax - yMin || 1;
    const zRange = zMax - zMin || 1;

    // Grid resolution
    const bins = is3D ? Math.ceil(Math.cbrt(threshLOD)) : Math.ceil(Math.sqrt(threshLOD));
    const gridSize = is3D ? bins * bins * bins : bins * bins;

    // Use Int32Array for performance
    const grid = new Int32Array(gridSize).fill(-1);

    // Pre-calculate factors
    const xFactor = bins / xRange;
    const yFactor = bins / yRange;
    const zFactor = bins / zRange;

    // Loop variables
    let xi: number, yi: number, zi: number, index: number;

    for (let j = 0; j < visibleIndices.length; j++) {
        const i = visibleIndices[j];
        const xVal = xValues[i];
        const yVal = yValues[i];

        // Map to grid
        xi = Math.floor((xVal - xMin) * xFactor);
        yi = Math.floor((yVal - yMin) * yFactor);

        // Clamp indices
        if (xi >= bins) xi = bins - 1;
        if (yi >= bins) yi = bins - 1;
        if (xi < 0) xi = 0;
        if (yi < 0) yi = 0;

        if (is3D && zValues) {
            zi = Math.floor((zValues[i] - zMin) * zFactor);
            if (zi >= bins) zi = bins - 1;
            if (zi < 0) zi = 0;

            // Flat index: x + y*w + z*w*h
            index = xi + yi * bins + zi * bins * bins;
        } else {
            index = xi + yi * bins;
        }

        if (grid[index] === -1) {
            grid[index] = i;
        }
    }

    const indices: number[] = [];
    for (let i = 0; i < grid.length; i++) {
        if (grid[i] !== -1) {
            indices.push(grid[i]);
        }
    }
    indices.sort((a, b) => a - b);

    return indices;
}
