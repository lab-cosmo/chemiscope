/**
 * @packageDocumentation
 * @module map
 */

import { Bounds, arrayMaxMin } from '../utils';
import { CameraState, projectPoints, getLookAtMatrix } from '../utils/camera';

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
 * @param priorityMask Optional boolean mask: ids with mask[id] === true win their
 *                     bin over ids with mask[id] === false (used to prefer
 *                     foreground points when a selection filter is active).
 * @returns Array of indices to display
 */
export function computeScreenSpaceLOD(
    xValues: number[],
    yValues: number[],
    zValues: number[],
    camera: CameraState,
    bounds?: Bounds,
    maxPoints: number = 50000,
    aspectRatio: number = 1.0,
    priorityMask?: boolean[]
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
        // Foreground ids (priorityMask[id] === true) get the high bit cleared,
        // background ids get it set, so any foreground id always beats any
        // background id competing for the same bin. The `>>> 0` keeps the
        // result an unsigned 32-bit value so it compares correctly against
        // bestHash, which is read from a Uint32Array.
        const baseH = Math.imul(id, 2654435761) >>> 1;
        const h = priorityMask && !priorityMask[id] ? (baseH | 0x80000000) >>> 0 : baseH;
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
 * @param priorityMask Optional boolean mask: ids with mask[id] === true win their
 *                     bin over ids with mask[id] === false (used to prefer
 *                     foreground points when a selection filter is active).
 * @returns Array of indices to display
 */
export function computeLODIndices(
    xValues: number[],
    yValues: number[],
    zValues: number[] | null,
    bounds?: Bounds,
    maxPoints: number = 50000,
    priorityMask?: boolean[]
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

        // Foreground ids (priorityMask[id] === true) get the high bit cleared,
        // background ids get it set, so any foreground id always beats any
        // background id competing for the same bin. The `>>> 0` keeps the
        // result an unsigned 32-bit value so it compares correctly against
        // bestHash, which is read from a Uint32Array.
        const baseH = Math.imul(id, 2654435761) >>> 1;
        const h = priorityMask && !priorityMask[id] ? (baseH | 0x80000000) >>> 0 : baseH;
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

// cap the finest grids at about one million cells in 2d and two million in 3d
const LEVELS: Record<number, number> = { 2: 10, 3: 7 };
// keep some points outside the current view
const GLOBAL_FRACTION = 0.2;
// cover small view changes while waiting for the next update
const WINDOW_PADDING = 0.15;

// axis limits: [[xMin, xMax], [yMin, yMax]], plus [zMin, zMax] in 3d
type Box = [number, number][];
type CameraPose = { x: number; y: number; z: number };

/**
 * Build a fixed order of points
 * Zoom and pan choose points from this order without rebuilding it
 * Only changes to axes or filters need a new sampler
 */
export class LODSampler {
    // ids of the points to display, sorted
    public indices: number[] = [];

    // number of axes
    private _dims: number;
    // minimum value on each axis across the full dataset
    private _min: number[] = [];
    // max - min on each axis or 1 if equal
    private _range: number[] = [];

    // per-axis values mapped to [0, 1], indexed by original point id
    // NaN marks a coordinate that cannot be plotted
    private _coords: Float64Array[] = [];

    // point ids in rank order, independent of the current view
    private _order: Uint32Array;
    private _maxPoints: number;
    private _globalCount: number;

    /**
     * @param coordinates x/y(/z) values of all points
     * @param maxPoints maximum number of points selected inside the view
     * @param priority points with "priority" flag win their cell over the others
     * @param visible points with "visible == false" are never selected
     */
    constructor(
        coordinates: number[][],
        maxPoints: number,
        priority?: boolean[],
        visible?: boolean[]
    ) {
        const n = coordinates[0].length;
        this._dims = coordinates.length;
        this._maxPoints = maxPoints;
        this._globalCount = Math.ceil(maxPoints * GLOBAL_FRACTION);

        // use each axis full-data min and max, not the current zoom limits
        // this keeps grid cells in the same places when the view changes
        for (const values of coordinates) {
            let min = Infinity;
            let max = -Infinity;

            for (const v of values) {
                if (isFinite(v)) {
                    min = Math.min(min, v);
                    max = Math.max(max, v);
                }
            }

            // if no value is finite, all coordinates below will be NaN
            if (!isFinite(min)) {
                min = 0;
                max = 1;
            }

            const range = max > min ? max - min : 1;
            const norm = new Float64Array(n);

            // map varying axes from min = 0 to max = 1 for the sampling grid
            // e.g. values 10, 20, 30 become 0, 0.5, 1
            for (let i = 0; i < n; i++) {
                norm[i] = isFinite(values[i]) ? (values[i] - min) / range : NaN;
            }

            this._min.push(min);
            this._range.push(range);
            this._coords.push(norm);
        }

        // only build the ranking once, later view changes reuse it
        this._order = this._rank(priority, visible);
    }

    /**
     * Select the view using the existing ranking, or the full range if omitted
     * Returns true if the selected points changed
     */
    public select(view?: Box): boolean {
        const window: Box = [];

        for (let d = 0; d < this._dims; d++) {
            let lo = 0;
            let hi = 1;

            if (view !== undefined) {
                // convert the axis limits to the [0, 1] coordinates in _coords
                lo = (Math.min(view[d][0], view[d][1]) - this._min[d]) / this._range[d];
                hi = (Math.max(view[d][0], view[d][1]) - this._min[d]) / this._range[d];

                // the displayed range can extend beyond the data
                lo = isNaN(lo) ? 0 : Math.min(1, Math.max(0, lo));
                hi = isNaN(hi) ? 1 : Math.min(1, Math.max(lo, hi));
            }

            // include nearby points so small pans do not expose empty edges
            const pad = (hi - lo) * WINDOW_PADDING;
            window.push([Math.max(0, lo - pad), Math.min(1, hi + pad)]);
        }

        // a point is inside the view only if every axis is within its padded limits
        return this._pick((id) => {
            for (let d = 0; d < this._dims; d++) {
                const u = this._coords[d][id];
                if (u < window[d][0] || u > window[d][1]) {
                    return false;
                }
            }
            return true;
        });
    }

    /**
     * Choose points inside the current 3D view using Plotly camera
     * Convert each point to its position in the scene, then to its screen position
     * Returns true if the chosen point ids changed
     *
     * @param bounds axis ranges of the scene
     * @param camera scene camera
     * @param aspect scene aspect ratio
     * @param widthOverHeight pixel aspect ratio of the scene
     */
    public selectCamera(
        bounds: Bounds,
        camera: { eye: CameraPose; center: CameraPose; up: CameraPose },
        aspect: CameraPose,
        widthOverHeight: number
    ): boolean {
        const view = getLookAtMatrix(camera);
        const ranges = [bounds.x, bounds.y, bounds.z as [number, number]];
        const scale = [aspect.x, aspect.y, aspect.z];

        // axis starts and widths
        const lo = ranges.map(([a], d) => (a - this._min[d]) / this._range[d]);

        // reversed axis has a negative width preserving its direction
        const span = ranges.map(([a, b], d) => (b - a) / this._range[d]);

        // the screen's top edge is at 1, add a margin for small view changes
        const limit = 1 + WINDOW_PADDING;

        // reuse this scratch space for each point
        const world = [0, 0, 0];

        // tell _pick whether this point would appear inside the padded screen
        return this._pick((id) => {
            for (let d = 0; d < 3; d++) {
                // position within this axis (0 at its start, 1 at its end)
                const t = (this._coords[d][id] - lo[d]) / span[d];

                // points outside the displayed axis limits are clipped by plotly
                if (!(t >= 0 && t <= 1)) {
                    return false;
                }

                // put the midpoint at 0, then stretch to plotly's axis length
                // with aspect 1, the axis runs from -0.5 to 0.5
                world[d] = (t - 0.5) * scale[d];
            }

            // apply the camera rotation and translation to get screen x and y
            const vx = view[0] * world[0] + view[4] * world[1] + view[8] * world[2] + view[12];
            const vy = view[1] * world[0] + view[5] * world[1] + view[9] * world[2] + view[13];

            // a wide plot shows more horizontally, both edges include the margin
            return Math.abs(vx) <= limit * widthOverHeight && Math.abs(vy) <= limit;
        });
    }

    /**
     * Walk the fixed order until the view's point budget is filled
     * Keep a few points outside the view so pans and rotations do not reveal empty space
     * Save the chosen ids and report whether the plot needs a redraw
     */
    private _pick(inWindow: (id: number) => boolean): boolean {
        const result: number[] = [];
        let inside = 0;

        // a smaller view keeps its existing points ahead of newly revealed detail
        for (let j = 0; j < this._order.length; j++) {
            const id = this._order[j];

            if (inWindow(id)) {
                result.push(id);

                // background points outside the view do not use this budget
                inside += 1;

                if (inside >= this._maxPoints) {
                    // every remaining point has a lower rank
                    break;
                }
            } else if (j < this._globalCount) {
                // keep a fixed background sample as the view moves
                result.push(id);
            }
        }

        // compare point sets in a consistent order to avoid needless redraws
        result.sort((a, b) => a - b);

        const changed =
            result.length !== this.indices.length || result.some((id, i) => id !== this.indices[i]);
        this.indices = result;
        return changed;
    }

    /**
     * Build the point order once using nested grids over the full data range
     * Each cell keeps the point with the smallest hash, as in the original sampler
     * Winners of larger cells come first, other points supply detail when zooming in
     */
    private _rank(priority?: boolean[], visible?: boolean[]): Uint32Array {
        const n = this._coords[0].length;
        const dims = this._dims;
        const levels = LEVELS[dims];

        // a hash avoids favoring points that happen to be first in the input
        // leave the high bit free for foreground priority
        const hash = new Uint32Array(n);
        for (let i = 0; i < n; i++) {
            const h = Math.imul(i, 2654435761) >>> 1;

            // setting the high bit makes background points lose within their cell
            hash[i] = priority && !priority[i] ? (h | 0x80000000) >>> 0 : h;
        }

        // a point needs valid coordinates on every axis to enter the ranking
        const eligible: number[] = [];
        for (let i = 0; i < n; i++) {
            if (visible && !visible[i]) {
                continue;
            }

            let valid = true;
            for (let d = 0; d < dims; d++) {
                if (isNaN(this._coords[d][i])) {
                    valid = false;
                    break;
                }
            }

            if (valid) {
                eligible.push(i);
            }
        }

        if (eligible.length === 0) {
            return new Uint32Array(0);
        }

        // lower levels rank first, points that never win start last
        // level zero is reserved for each axis's minimum and maximum
        const level = new Uint8Array(n).fill(levels + 2);
        const cells = 1 << (levels * dims);

        // allocate for the finest grid and reuse its storage for coarser ones
        const bestHash = new Uint32Array(cells);
        const bestId = new Int32Array(cells);
        let candidates = eligible;

        // start with small cells, then halve the number of cells along each axis
        // each larger cell chooses a winner from its smaller cells' winners
        for (let k = levels; k >= 0; k--) {
            const size = 1 << k;
            const count = size ** dims;

            // clear just the cells used by this level
            bestHash.fill(0xffffffff, 0, count);
            bestId.fill(-1, 0, count);

            for (const id of candidates) {
                let idx = 0;

                // flatten the cell coordinates into one array index
                // a point at 1 belongs to the last cell
                for (let d = dims - 1; d >= 0; d--) {
                    idx = idx * size + Math.min(size - 1, Math.floor(this._coords[d][id] * size));
                }
                if (hash[id] < bestHash[idx]) {
                    bestHash[idx] = hash[id];
                    bestId[idx] = id;
                }
            }

            // only cell winners can advance to the next, coarser grid
            const winners: number[] = [];
            for (let c = 0; c < count; c++) {
                if (bestId[c] >= 0) {
                    winners.push(bestId[c]);
                    // a coarser win moves this point earlier in the ranking
                    level[bestId[c]] = k + 1;
                }
            }

            // losing points stay in eligible for later detail
            candidates = winners;
        }

        // plotly autoscales from the sampled points, not the full dataset
        for (let d = 0; d < dims; d++) {
            let minIndex = eligible[0];
            let maxIndex = eligible[0];

            // find the point ids with the smallest and largest value on this axis
            for (const id of eligible) {
                const u = this._coords[d][id];
                if (u < this._coords[d][minIndex]) {
                    minIndex = id;
                }
                if (u > this._coords[d][maxIndex]) {
                    maxIndex = id;
                }
            }

            // rank both endpoints before all grid winners
            level[minIndex] = 0;
            level[maxIndex] = 0;
        }

        // coarse-grid winners come first so sparse regions stay represented
        // within each level, prefer foreground points
        eligible.sort((a, b) => level[a] - level[b] || hash[a] - hash[b] || a - b);
        return Uint32Array.from(eligible);
    }
}
