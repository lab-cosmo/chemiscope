/**
 * @packageDocumentation
 * @module map
 */

import { arrayMaxMin } from '../utils';
import { CameraState, getLookAtMatrix, projectPoints } from '../utils/camera';

/**
 * Computes the subset of points to display based on 2D screen-space projection.
 * This is used for 3D plots where the visible density depends on camera rotation/zoom.
 *
 * @param xValues Array of X coordinates
 * @param yValues Array of Y coordinates
 * @param zValues Array of Z coordinates
 * @param camera  Current camera state (eye, center, up, zoom)
 * @param threshLOD Maximum number of points to display (default: 50000)
 * @returns Sorted array of indices to display
 */
export function computeScreenSpaceLOD(
    xValues: number[],
    yValues: number[],
    zValues: number[],
    camera: CameraState,
    bounds: { x: [number, number]; y: [number, number]; z?: [number, number] },
    threshLOD: number = 50000
): number[] {
    const nPoints = xValues.length;
    // Calculate View Matrix once
    const viewMatrix = getLookAtMatrix(camera);
    const projections = projectPoints(xValues, yValues, zValues, camera, bounds);

    // Grid resolution: we want roughly threshLOD points on screen.
    const bins = Math.ceil(Math.sqrt(threshLOD));
    const grid = new Map<string, number>();

    // this is the range we use to bin points
    const viewSize = 2.0;    
    const uStep = viewSize / bins;
    const vStep = viewSize / bins;

    // Use a slightly larger clip bound to avoid popping at edges
    const clipSize = viewSize * 1.1;

    for (let i = 0; i < nPoints; i++) {
        const u = projections.x[i];
        const v = projections.y[i];

        // Clip points well outside the view
        if (Math.abs(u) > clipSize || Math.abs(v) > clipSize) {
            continue;
        }

        // Center the grid around 0
        const ui = Math.floor(u / uStep);
        const vi = Math.floor(v / vStep);

        const key = `${ui}_${vi}`;
        if (!grid.has(key)) {
            grid.set(key, i);
        }
    }

    const indices = Array.from(grid.values());
    indices.sort((a, b) => a - b);
    console.log("selected pooints: ", indices.length)
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
            // Default Z range if not provided in 3D bounds (unlikely but safe)
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

    // Avoid division by zero
    const xRange = xMax - xMin || 1;
    const yRange = yMax - yMin || 1;
    const zRange = zMax - zMin || 1;

    // Grid resolution, determined so that for a dense dataset we get roughly threshLOD points
    const bins = is3D ? Math.ceil(Math.cbrt(threshLOD)) : Math.ceil(Math.sqrt(threshLOD));
    const grid = new Map<string, number>();
    const clip = Array<number>(0);

    // Re-use loop variables for performance
    let xi: number, yi: number, zi: number, key: string;

    for (let i = 0; i < nPoints; i++) {
        const xVal = xValues[i];
        const yVal = yValues[i];

        // Clipping: strictly ignore points outside the view if bounds are provided
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
        clip.push(i); // Keep track of points inside the bounds

        // Calculate grid coordinates relative to the current View/Range
        xi = Math.floor(((xVal - xMin) / xRange) * bins);
        yi = Math.floor(((yVal - yMin) / yRange) * bins);

        if (is3D && zValues) {
            zi = Math.floor(((zValues[i] - zMin) / zRange) * bins);
            key = `${xi}_${yi}_${zi}`;
        } else {
            key = `${xi}_${yi}`;
        }

        // Only store the first point found in this grid cell
        if (!grid.has(key)) {
            grid.set(key, i);
        }
    }

    if (clip.length < threshLOD) {
        // If the number of points inside the bounds is already below the threshold,
        // return all those points without further downsampling.
        return clip;
    }

    // Convert map to sorted array of indices
    const indices = Array.from(grid.values());
    indices.sort((a, b) => a - b);

    return indices;
}
