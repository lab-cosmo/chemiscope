/**
 * @packageDocumentation
 * @module map
 */

import { arrayMaxMin } from '../utils';
import { CameraState, projectPoints } from '../utils/camera';
import { NumericProperty } from './data';
import { MapOptions } from './options';

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
function computeScreenSpaceLOD(
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

    // Grid resolution: we want roughly threshLOD points on screen.
    const bins = Math.ceil(Math.sqrt(threshLOD));
    const grid = new Map<string, number>();

    // this is the range we use to bin points
    const viewSize = 2.0;
    // Use a slightly larger clip bound to avoid popping at edges
    const clipSize = viewSize * 1.1;
    const uStep = (clipSize * 2) / bins;
    const vStep = (clipSize * 2) / bins;

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
function computeLODIndices(
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

interface LODDependencies {
    xProperty: string;
    yProperty: string;
    zProperty: string;
    lodEnabled: boolean;
}

/**
 * High level controller managing LOD indices, debouncing and restyling.
 */
export class LODManager {
    private _indices: number[] | null = null;
    private _lock = false;
    private _debounceTimer: number | undefined;
    private _lastDependencies: LODDependencies | null = null;
    private _cachedValues: {
        x: number[];
        y: number[];
        z: number[] | null;
    } | null = null;

    constructor(
        private options: MapOptions,
        private is3D: () => boolean,
        private getProperty: (name: string) => NumericProperty,
        private restyleFull: () => Promise<void>,
        private threshLOD: number = 50000,
        private debounceMs: number = 300,
        private settleMs: number = 200
    ) {}

    public get indices(): number[] | null {
        return this._indices;
    }

    public isLocked(): boolean {
        return this._lock;
    }

    public setLock(v: boolean): void {
        this._lock = v;
    }

    private _applyLODNumeric(values: number[]): number[] {
        if (this._indices === null) {
            return values;
        }

        const len = this._indices.length;
        const result = new Float64Array(len);
        for (let i = 0; i < len; i++) {
            result[i] = values[this._indices[i]];
        }

        return Array.from(result);
    }

    public applyLOD<T>(values: T[]): T[] {
        if (this._indices === null) {
            return values;
        }

        if (values.length > 0 && typeof values[0] === 'number') {
            return this._applyLODNumeric(values as unknown as number[]) as unknown as T[];
        }

        const len = this._indices.length;
        const result = new Array<T>(len);
        for (let i = 0; i < len; i++) {
            result[i] = values[this._indices[i]];
        }

        return result;
    }

    private _getCurrentDependencies(): LODDependencies {
        return {
            xProperty: this.options.x.property.value,
            yProperty: this.options.y.property.value,
            zProperty: this.options.z.property.value,
            lodEnabled: this.options.useLOD.value,
        };
    }

    private _needsRecomputation(): boolean {
        const current = this._getCurrentDependencies();

        if (this._lastDependencies === null) {
            return true; // First time
        }

        return (
            current.xProperty !== this._lastDependencies.xProperty ||
            current.yProperty !== this._lastDependencies.yProperty ||
            current.zProperty !== this._lastDependencies.zProperty ||
            current.lodEnabled !== this._lastDependencies.lodEnabled
        );
    }

    public computeLOD(bounds?: {
        x: [number, number];
        y: [number, number];
        z?: [number, number];
    }): void {
        if (!this._needsRecomputation() && bounds === undefined) {
            return;
        }

        this._lastDependencies = this._getCurrentDependencies();

        if (!this.options.useLOD.value) {
            this._indices = null;
            this._cachedValues = null;
            return;
        }

        const xProp = this.options.x.property.value;
        const yProp = this.options.y.property.value;
        const zProp = this.options.z.property.value;

        const xValues = this.getProperty(xProp).values;

        if (xValues.length <= this.threshLOD) {
            this._indices = null;
            this._cachedValues = null;
            return;
        }

        const yValues = this.getProperty(yProp).values;
        const is3D = this.is3D() && zProp !== '';
        const zValues = is3D ? this.getProperty(zProp).values : null;

        this._cachedValues = {
            x: xValues,
            y: yValues,
            z: zValues,
        };

        const lodIndices = computeLODIndices(
            xValues,
            yValues,
            zValues,
            undefined,
            Math.floor(this.threshLOD / 10)
        );

        if (is3D && zValues && this.options.camera.value && bounds) {
            lodIndices.push(
                ...computeScreenSpaceLOD(
                    xValues,
                    yValues,
                    zValues,
                    this.options.camera.value,
                    bounds,
                    this.threshLOD
                )
            );
        } else {
            lodIndices.push(
                ...computeLODIndices(xValues, yValues, zValues, bounds, this.threshLOD)
            );
        }

        this._indices = [...new Set(lodIndices)].sort((a, b) => a - b);
    }

    public forceComputeLOD(bounds?: {
        x: [number, number];
        y: [number, number];
        z?: [number, number];
    }): void {
        this._lastDependencies = null; // Force recomputation
        this.computeLOD(bounds);
    }

    public async restyleLOD(): Promise<void> {
        await this.restyleFull();
    }

    public scheduleLODUpdate(bounds: {
        x: [number, number];
        y: [number, number];
        z?: [number, number];
    }): void {
        if (this._debounceTimer !== undefined) {
            window.clearTimeout(this._debounceTimer);
        }

        this._debounceTimer = window.setTimeout(() => {
            this._debounceTimer = undefined;
            void this.performLODUpdate(bounds);
        }, this.debounceMs);
    }

    private async performLODUpdate(bounds: {
        x: [number, number];
        y: [number, number];
        z?: [number, number];
    }): Promise<void> {
        if (this._lock) {
            return;
        }

        this._lock = true;

        try {
            this.computeLOD(bounds);
            await this.restyleFull();
            await new Promise((resolve) => setTimeout(resolve, this.settleMs));
        } catch (error) {
            // avoid breaking the UI loop
            // eslint-disable-next-line no-console
            console.error('LOD update failed:', error);
        } finally {
            this._lock = false;
        }
    }

    public cancelPending(): void {
        if (this._debounceTimer !== undefined) {
            window.clearTimeout(this._debounceTimer);
            this._debounceTimer = undefined;
        }
    }

    public clearCache(): void {
        this._cachedValues = null;
    }
}
