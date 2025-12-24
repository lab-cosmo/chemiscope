/**
 * @packageDocumentation
 * @module map
 */

import assert from 'assert';

import Plotly from './plotly/plotly-scatter';
import { Config, Data, Layout, PlotlyScatterElement } from './plotly/plotly-scatter';
import * as plotlyStyles from './plotly/plotly-styles';
import fixPlot from './plotly/fix-plot';

import { Property, Settings } from '../dataset';

import { DisplayTarget, EnvironmentIndexer, Indexes } from '../indexer';
import { OptionModificationOrigin } from '../options';
import { GUID, PositioningCallback, Warnings, arrayMaxMin } from '../utils';
import { enumerate, getElement, getFirstKey } from '../utils';

import { MapData, NumericProperties, NumericProperty } from './data';
import { MarkerData } from './marker';
import { AxisOptions, MapOptions, get3DSymbol } from './options';
import * as styles from '../styles';

import PNG_SVG from '../static/download-png.svg';
import SVG_SVG from '../static/download-svg.svg';

const DEFAULT_LAYOUT = {
    // coloraxis is used for the markers
    coloraxis: {
        cmax: 0,
        cmin: 0,
        colorbar: {
            len: 1,
            thickness: 20,
            title: {
                text: '',
                side: 'right',
                font: {
                    size: 15,
                },
            },
            y: 0,
            yanchor: 'bottom',
        },
        colorscale: [] as Plotly.ColorScale,
        showscale: true,
    },
    hovermode: 'closest',
    legend: {
        itemclick: false,
        itemdoubleclick: false,
        tracegroupgap: 5,
        y: 1,
        yanchor: 'top',
    },
    margin: {
        b: 50,
        l: 50,
        r: 50,
        t: 50,
    },
    scene: {
        camera: {
            projection: {
                type: 'orthographic',
            },
        },
        xaxis: {
            showspikes: false,
            title: '',
        },
        yaxis: {
            showspikes: false,
            title: '',
        },
        zaxis: {
            showspikes: false,
            title: '' as undefined | string,
        },
    },
    showlegend: true,
    xaxis: {
        range: undefined as (number | undefined)[] | undefined,
        title: '',
        type: 'linear',
        zeroline: false,
    },
    yaxis: {
        range: undefined as (number | undefined)[] | undefined,
        title: '',
        type: 'linear',
        zeroline: false,
    },
    zaxis: {
        range: undefined as (number | undefined)[] | undefined,
        title: '',
        type: 'linear',
        zeroline: false,
    },
};

const DEFAULT_CONFIG = {
    displayModeBar: true,
    displaylogo: false,
    responsive: true,
    scrollZoom: true,

    modeBarButtonsToRemove: [
        'hoverClosestCartesian',
        'hoverCompareCartesian',
        'toggleSpikelines',
        'autoScale2d',
        'zoomIn2d',
        'zoomOut2d',
        'select2d',
        'lasso2d',
        'hoverClosest3d',
        'tableRotation',
        'resetCameraLastSave3d',
        'toImage',
    ],

    modeBarButtonsToAdd: [
        [
            {
                name: 'Download PNG',
                icon: {
                    width: 400,
                    height: 447,
                    path: extractSvgPath(PNG_SVG),
                },
                click: function (gd: PlotlyScatterElement) {
                    const width = Math.max(gd._fullLayout.width, 600);
                    const ratio = gd._fullLayout.height / gd._fullLayout.width;
                    const height = width * ratio;

                    Plotly.downloadImage(gd, {
                        filename: 'chemiscope-map',
                        format: 'png',
                        width: width,
                        height: height,
                        // scale is not part of `DownloadImgopts`, but accepted
                        // by the function anyway
                        scale: 3,
                    } as unknown as Plotly.DownloadImgopts).catch((e: unknown) =>
                        setTimeout(() => {
                            throw e;
                        })
                    );
                },
            },
        ],
        [
            {
                name: 'Download SVG',
                icon: {
                    width: 400,
                    height: 447,
                    path: extractSvgPath(SVG_SVG),
                },
                click: function (gd: PlotlyScatterElement) {
                    Plotly.downloadImage(gd, {
                        filename: 'chemiscope-map',
                        format: 'svg',
                        width: Math.max(gd._fullLayout.width, 600),
                        height: Math.max(gd._fullLayout.height, 600),
                    }).catch((e: unknown) =>
                        setTimeout(() => {
                            throw e;
                        })
                    );
                },
            },
        ],
    ],
};

/**
 * The {@link PropertiesMap} class displays a 2D or 3D map (scatter plot) of
 * properties in the dataset, using [plotly.js](https://plot.ly/javascript/)
 * for rendering.
 *
 * Properties can be used as x, y, or z values, as well as points color and
 * size. Additionally, string properties can be used as symbols for the scatter
 * plot markers.
 */
export class PropertiesMap {
    /** Callback fired when the plot is clicked and the position of the active marker changes */
    public onselect: (indexes: Indexes) => void;
    /**
     * Callback fired when the active marker is changed by clicking on the map
     *
     * @param guid GUID of the new active marker
     * @param indexes index of the environment the new active marker is showing
     */
    public activeChanged: (guid: GUID, indexes: Indexes) => void;

    /**
     * Callback to get the initial positioning of the settings modal.
     *
     * The callback gets the current placement of the settings as a
     * [DOMRect](https://developer.mozilla.org/en-US/docs/Web/API/DOMRect), and
     * should return top and left positions in pixels, used with `position:
     * fixed`. The callback is called once, the first time the settings are
     * opened.
     */
    public positionSettingsModal: PositioningCallback;

    public warnings: Warnings;

    /// Shadow root for isolation
    private _shadow: ShadowRoot;
    /// HTML element holding the full plot
    private _root: HTMLElement;
    /// Plotly plot
    private _plot!: PlotlyScatterElement;
    /// All known properties
    private _data: MapData;

    /// GUID of the currently selected point
    private _active?: GUID;
    /// Map associating currently selected markers GUID to additional data
    private _selected: Map<GUID, MarkerData>;

    /// environment indexer
    private _indexer: EnvironmentIndexer;
    // widget display target
    private _target: DisplayTarget;
    /// Settings of the map
    private _options: MapOptions;
    /// Button used to reset the range of color axis
    private _colorReset: HTMLButtonElement;
    /// Plotly fix instance
    private _plotFix!: ReturnType<typeof fixPlot>;

    private _lodEnabled = true;
    private _lodMaxPoints = 30000; // tune 20000–50000
    private _lodMap: number[] = []; // display index -> environment index
    private _lodDirty = true;
    private _lodTimer: number | undefined;

    private _lodCache:
    | {
          x: number[];
          y: number[];
          color: Array<string | number>;
          size: Array<number>;
          symbol: string[] | number[] | string; // depends on your symbol mode
      }
    | undefined;

    private _pickScreenBinnedIndices(
        x: number[],
        y: number[],
        maxPoints: number
    ): number[] {
        // Plotly internal axis objects (works in your code already via _pixelCoordinate usage)
        const xaxis = (this._plot as any)._fullLayout?.xaxis;
        const yaxis = (this._plot as any)._fullLayout?.yaxis;
        if (!xaxis || !yaxis) {
            // fallback: first maxPoints points
            const n = Math.min(x.length, maxPoints);
            return Array.from({ length: n }, (_, i) => i);
        }

        const w: number = xaxis._length; // pixels
        const h: number = yaxis._length; // pixels
        if (!w || !h) return [];

        // Choose a pixel bin size so total bins ≈ maxPoints
        const binPx = Math.max(2, Math.round(Math.sqrt((w * h) / maxPoints)));
        const binsX = Math.max(1, Math.floor(w / binPx));
        const binsY = Math.max(1, Math.floor(h / binPx));
        const binsN = binsX * binsY;

        // store chosen point per bin, and optionally choose the point closest to bin center
        const chosen = new Int32Array(binsN);
        chosen.fill(-1);
        const best = new Float32Array(binsN);
        best.fill(Number.POSITIVE_INFINITY);

        for (let i = 0; i < x.length; i++) {
            const xi = x[i];
            const yi = y[i];

            // l2p returns pixel coordinate in plot area coordinates for linear axes
            // (this is what your code already relies on)
            const px = xaxis.l2p(xi);
            const py = yaxis.l2p(yi);
            if (!Number.isFinite(px) || !Number.isFinite(py)) continue;
            if (px < 0 || px > w || py < 0 || py > h) continue;

            const bx = Math.min(binsX - 1, Math.max(0, (px / binPx) | 0));
            const by = Math.min(binsY - 1, Math.max(0, (py / binPx) | 0));
            const b = bx + by * binsX;

            // keep the point closest to the bin center for nicer visuals
            const cx = (bx + 0.5) * binPx;
            const cy = (by + 0.5) * binPx;
            const d2 = (px - cx) * (px - cx) + (py - cy) * (py - cy);
            if (d2 < best[b]) {
                best[b] = d2;
                chosen[b] = i;
            }
        }

        const out: number[] = [];
        for (let b = 0; b < binsN; b++) {
            const i = chosen[b];
            if (i >= 0) out.push(i);
        }
        return out;
    }

    private _sliceNumberArray(full: number[], idx: number[]): number[] {
        const out = new Array<number>(idx.length);
        for (let i = 0; i < idx.length; i++) out[i] = full[idx[i]];
        return out;
    }

    private _sliceColorArray(
        full: Array<string | number>,
        idx: number[]
    ): Array<string | number> {
        const out = new Array<string | number>(idx.length);
        for (let i = 0; i < idx.length; i++) out[i] = full[idx[i]];
        return out;
    }

    private _sliceSymbolArray(
        full: string[] | number[] | string,
        idx: number[]
    ): string[] | number[] | string {
        // If symbol is a scalar (e.g. 'circle'), keep it scalar
        if (typeof full === 'string') return full;

        const out = new Array<any>(idx.length);
        for (let i = 0; i < idx.length; i++) out[i] = (full as any)[idx[i]];
        return out as any;
    }

    private _rebuildLodCache(): void {
        // FULL arrays for trace 0 (main)
        const x = this._coordinates(this._options.x, 0)[0] as number[];
        const y = this._coordinates(this._options.y, 0)[0] as number[];

        const color = this._colors(0)[0]; // full marker.color array
        const size = this._sizes(0)[0] as number[]; // full marker.size array
        const symbol = this._symbols(0)[0]; // can be scalar or array

        this._lodCache = { x, y, color, size, symbol };
        this._lodDirty = false;
    }    

    private _scheduleLodUpdate(): void {
        if (!this._lodEnabled) return;
        if (this._options.joinPoints.value) return; // LOD + lines is usually wrong

        if (this._lodTimer !== undefined) {
            window.clearTimeout(this._lodTimer);
        }
        this._lodTimer = window.setTimeout(() => this._applyLodFromView(), 0);
    }

    private _applyLodFromView(): void {
        if (this._is3D()) return;
        if (!this._lodEnabled) return;
        if (this._options.joinPoints.value) return;

        if (this._lodDirty || !this._lodCache) {
            this._rebuildLodCache();
        }
        const cache = this._lodCache!;
        const n = cache.x.length;

        // if small enough, just show everything (and clear mapping)
        if (n <= this._lodMaxPoints) {
            this._lodMap = [];
            // ensure trace 0 is full if we previously downsampled
            Plotly.restyle(
                this._plot,
                {
                    x: [cache.x],
                    y: [cache.y],
                    'marker.color': [cache.color],
                    'marker.size': [cache.size],
                    'marker.symbol': [cache.symbol],
                },
                [0]
            ).catch(() => {});
            return;
        }

        const idx = this._pickScreenBinnedIndices(cache.x, cache.y, this._lodMaxPoints);
        this._lodMap = idx; // display->original

        const xs = this._sliceNumberArray(cache.x, idx);
        const ys = this._sliceNumberArray(cache.y, idx);
        const cs = this._sliceColorArray(cache.color, idx);
        const ss = this._sliceNumberArray(cache.size, idx);
        const sy = this._sliceSymbolArray(cache.symbol, idx);

        Plotly.restyle(
            this._plot,
            {
                x: [xs],
                y: [ys],
                'marker.color': [cs],
                'marker.size': [ss],
                'marker.symbol': [sy],
            },
            [0]
        ).catch(() => {});
    }


    /**
     * Create a new {@link PropertiesMap} inside the DOM element with the given HTML
     * `id`
     *
     * @param element    HTML element or string 'id' of the element where
     *                   the map should live
     * @param settings   settings for all panels
     * @param indexer    {@link EnvironmentIndexer} used to translate indexes from
     *                   environments index to structure/atom indexes
     * @param target       widget display target, either stucture or atom
     * @param properties properties to be displayed
     */
    constructor(
        element: string | HTMLElement,
        settings: Settings,
        indexer: EnvironmentIndexer,
        target: DisplayTarget,
        properties: { [name: string]: Property },
        warnings: Warnings
    ) {
        this._indexer = indexer;
        this._target = target;
        this.onselect = () => {};
        this.activeChanged = () => {};
        this._selected = new Map<GUID, MarkerData>();

        this.warnings = warnings;

        // DOM structure outside the map:
        // - containerElement/element (#chemiscope-map)
        //   - hostElement (layer needed for removal)
        //     - this._shadow (shadow root)
        //       - this._root
        //         - plotly root
        //         - viewer button
        //         - marker

        // Attach a shadow DOM to the host element for isolation
        const containerElement = getElement(element);
        const hostElement = document.createElement('div');
        hostElement.style.setProperty('height', '100%');
        containerElement.appendChild(hostElement);
        this._shadow = hostElement.attachShadow({ mode: 'open' });

        // Create a root element inside the shadow DOM to append plot to it
        this._root = document.createElement('div');
        this._root.style.setProperty('height', '100%');
        this._shadow.appendChild(this._root);
        if (this._root.style.position === '') {
            this._root.style.position = 'relative';
        }

        // Create a div for the plot
        this._plot = document.createElement('div') as unknown as PlotlyScatterElement;
        this._plot.style.width = '100%';
        this._plot.style.height = '100%';
        this._root.appendChild(this._plot);

        // Initialize data with the given properties
        this._data = new MapData(properties, this.warnings);

        // Initialize options used in the modal
        const currentProperties = this._getCurrentProperties();
        this._options = new MapOptions(
            this._root,
            currentProperties,
            (rect) => this.positionSettingsModal(rect),
            settings,
            this.warnings
        );
        this._colorReset = this._options.getModalElement<HTMLButtonElement>('map-color-reset');

        // Connect the settings to event listeners or handlers
        this._connectSettings();

        // Define the default position for the settings modal on top of the plot,
        // centered horizontally
        this.positionSettingsModal = (rect: DOMRect) => {
            const rootRect = this._root.getBoundingClientRect();
            return {
                left: rootRect.left + rootRect.width / 2 - rect.width / 2,
                top: rootRect.top + 20,
            };
        };

        // Create the Plotly plot within the plot element
        this._createPlot();

        // Adopt styles with the plot stylesheets as last one because the plot
        // needs to be created to obtain it
        this._shadow.adoptedStyleSheets = [
            styles.bootstrap,
            styles.chemiscope,
            plotlyStyles.globalStyleSheet,
            plotlyStyles.getPlotStyleSheet(this._plot),
        ];
    }

    /**
     * Change display target and adapt the element to the new target
     * @param target display target
     */
    public async switchTarget(target: DisplayTarget): Promise<void> {
        // Check if the target value actually changed
        if (target !== this._target) {
            // Set new widget target
            this._target = target;

            // Process markers to correctly set then to the updated map
            this._handleMarkers();
            if (this._active !== undefined) {
                const activeMarker = this._selected.get(this._active);
                if (activeMarker !== undefined) {
                    this.onselect(
                        this._indexer.fromEnvironment(activeMarker.current, this._target)
                    );
                }
            }

            // Update the map options based on the chosen target
            this._setupMapOptions();

            // Append the callbacks to the new options
            this._connectSettings();

            // Re-render the plot with the new data and layout
            await this._react(this._getTraces(), this._getLayout());
        }
    }

    /**
     * Processes the markers by updating their current value to structure or environment index on the current target
     * In 'structure' target, it set structure index, and in 'atom' target its environment index
     */
    private _handleMarkers() {
        for (const [key, marker] of this._selected.entries()) {
            if (this._target === 'structure') {
                marker.current = this._indexer.fromEnvironment(marker.current, 'atom').structure;
            } else {
                const environment = this._indexer.fromStructure(
                    marker.current,
                    'atom'
                )?.environment;
                assert(environment !== undefined);
                marker.current = environment;
            }
            this._selected.set(key, marker);
        }
    }

    /**
     * Remove all HTML added by this {@link PropertiesMap} in the current document
     */
    public remove(): void {
        // Remove the the shadow root's host. It is not possible to remove the shadow root directly.
        this._shadow.host.remove();

        // Remove options
        this._options.remove();

        // remove SVG element created by Plotly
        document.getElementById('js-plotly-tester')?.remove();

        // Remove listeners on the document caused by the plot fix
        this._plotFix.disable();
    }

    /**
     * Change the environment indicated by the currently active marker to
     * the one with the given `indexes`
     */
    public select(indexes: Indexes): void {
        if (this._active === undefined) {
            throw Error('tries to update selected environment, but there is no active marker');
        }

        // data.select needs `indexes.environment`, so make sure it is defined
        if (indexes.environment === undefined) {
            const fullIndexes = this._indexer.fromStructureAtom(
                this._target,
                indexes.structure,
                indexes.atom
            );
            if (fullIndexes === undefined) {
                const atomStr = indexes.atom === undefined ? '' : ` / atom ${indexes.atom}`;
                throw Error(
                    `can not find the environnement for structure ${indexes.structure}` + atomStr
                );
            }
            indexes = fullIndexes;
        }

        const data = this._selected.get(this._active);
        assert(data !== undefined);
        // this prevents infinite recursion
        // will return false if the data already corresponds to this index
        if (data.select(indexes)) {
            this._updateMarkers();
        }
    }

    /**
     * Set the marker with given GUID as the active marker.
     *
     * @param guid the GUID of the new active viewer
     */
    public setActive(guid: GUID): void {
        if (this._active !== undefined) {
            const oldData = this._selected.get(this._active);
            assert(oldData !== undefined);
            oldData.deactivate();
        }

        this._active = guid;
        const data = this._selected.get(this._active);
        assert(data !== undefined);
        data.activate();

        if (this._is3D()) {
            this._restyle({ 'marker.size': this._sizes(1) } as Data, 1);
        }
    }

    /**
     * Add a new marker to the map. The new marker is set as the active one
     *
     * @param guid GUID of the new marker
     * @param color GUID of the marker indicating the new viewer
     * @param indexes indexes of the environment that the new marker should show
     */
    public addMarker(guid: GUID, color: string, indexes: Indexes): void {
        assert(!this._selected.has(guid));

        const data = new MarkerData(guid, color, indexes.environment);
        this._root.appendChild(data.marker);
        data.marker.onclick = () => {
            this.setActive(guid);
            this.activeChanged(guid, this._indexer.fromEnvironment(data.current, this._target));
        };
        this._selected.set(guid, data);
        this._updateMarkers([data]);
        this.setActive(guid);
    }

    /**
     * Removes a marker from the map.
     *
     * @param guid GUID of the marker to remove
     */
    public removeMarker(guid: GUID): void {
        // If we remove the active marker, let's change which one is active
        if (this._active === guid) {
            if (this._selected.size === 1) {
                // we are removing the last one
                this._active = undefined;
            } else {
                this.setActive(getFirstKey(this._selected, guid));
            }
        }

        // remove HTML marker
        const data = this._selected.get(guid);
        assert(data !== undefined);
        data.remove();

        this._selected.delete(guid);
        this._updateMarkers();
    }

    /**
     * Apply saved settings to the map.
     */
    public applySettings(settings: Settings): void {
        this._options.applySettings(settings);
    }

    /**
     * Save the values of the current settings in a way that an be used with
     * {@link applySettings} or saved to JSON.
     */
    public saveSettings(): Settings {
        return this._options.saveSettings();
    }

    /**
     * Add the given `callback` to be called whenever a setting changes. The
     * callback will be given the path to the settings as a list of keys; and
     * the new value of the setting.
     *
     * There is currently no way to remove a callback.
     */
    public onSettingChange(callback: (keys: string[], value: unknown) => void): void {
        this._options.onSettingChange(callback);
    }

    /**
     * Build the traces from the options data
     */
    private _getTraces(): Plotly.Data[] {
        const type = this._is3D() ? 'scatter3d' : 'scattergl';
        // The main trace, containing default data
        const main = {
            name: '',
            type: type,

            x: this._coordinates(this._options.x, 0)[0],
            y: this._coordinates(this._options.y, 0)[0],
            z: this._coordinates(this._options.z, 0)[0],

            hovertemplate: this._options.hovertemplate(),
            marker: {
                color: this._colors(0)[0],
                coloraxis: 'coloraxis',
                line: {
                    color: 'black',
                    width: this._options.markerOutline.value ? 0.5 : 0,
                },
                // prevent plolty from messing with opacity when doing bubble
                // style charts (different sizes for each point)
                opacity: this._options.color.opacity.value / 100,
                size: this._sizes(0)[0],
                sizemode: 'area',
                symbol: this._symbols(0)[0],
            },
            line: {
                // lines type (if required)
                color: 'black',
                width: 0.5,
                dash: 'solid',
            },
            mode: this._options.joinPoints.value ? 'lines+markers' : 'markers',
            showlegend: false,
        };

        // Create a second trace to store the last clicked point, in order to
        // display it on top of the main plot with different styling. This is
        // only used in 3D mode, since it is way slower than moving
        // this._selectedMarker around.
        const selected = {
            name: 'selected',
            type: type,

            x: [],
            y: [],
            z: [],

            hoverinfo: 'none',
            marker: {
                color: [],
                line: {
                    color: [],
                    width: 2,
                },
                opacity: 1,
                size: [],
                sizemode: 'area',
            },
            mode: 'markers',
            showlegend: false,
        };

        const traces = [main as Data, selected as Data];

        // Calculate legend names and show legend flags based on data properties
        const legendNames = this._legendNames().slice(2);
        const showlegend = this._showlegend().slice(2);
        assert(legendNames.length === showlegend.length);
        const currentLength = legendNames.length;

        // Adjust arrays based on maxSymbols configuration
        if (this._data.maxSymbols > 0) {
            // resize & fill arrays
            legendNames.length = this._data.maxSymbols;
            legendNames.fill('', currentLength);
            showlegend.length = this._data.maxSymbols;
            showlegend.fill(false, currentLength);
        }

        // add empty traces to be able to display the symbols legend
        // one trace for each possible symbol
        for (let s = 0; s < this._data.maxSymbols; s++) {
            const data = {
                name: legendNames[s],
                type: type,

                // We need to add a dummy point to force plotly to display the
                // associated legend; but we don't want to see the point in the
                // map. Setting the coordinates to NaN achieve this.
                x: [NaN],
                y: [NaN],
                z: [NaN],

                marker: {
                    color: 'black',
                    size: 10,
                    symbol: this._is3D() ? get3DSymbol(s) : s,
                },
                mode: 'markers',
                showlegend: showlegend[s],
            };
            traces.push(data as Data);
        }
        return traces;
    }

    /**
     * Update options with the structure or atom default or config settings
     */
    private _setupMapOptions() {
        // Helper function to update `_option` values with related target configuration
        const refreshMapOptions = (properties: NumericProperties) => {
            // Check if target's properties actually exist
            if (properties && Object.keys(properties).length > 1) {
                // Delete previous target modal from DOM
                this._options.remove();

                // Re-create modal
                this._options = new MapOptions(
                    this._root,
                    this._getCurrentProperties(),
                    (rect) => this.positionSettingsModal(rect),
                    {},
                    this.warnings
                );
            }
        };

        // Apply structure or atom options based on target
        if (this._target !== 'atom') {
            refreshMapOptions(this._data['structure']);
        } else {
            refreshMapOptions(this._data['atom']);
        }
    }

    /** Returns the properties related to the target (structure or environment) */
    private _getCurrentProperties() {
        const properties = this._data[this._target];
        const propertiesNames = Object.keys(properties);
        if (propertiesNames.length < 2) {
            // better error message in case the user forgot to give the
            // environments in the data
            if (this._target === 'structure') {
                if (Object.keys(this._data['atom']).length >= 2) {
                    throw Error(
                        'could not find enough structure properties to display, \
                        but there are atom properties. You could choose "atom" as \
                        target to display them.'
                    );
                }
            }

            let message = 'we need at least two properties to plot in the map';
            if (propertiesNames.length === 0) {
                message += ', we have none';
            } else {
                message += `, we have only one: '${propertiesNames[0]}'`;
            }

            throw Error(message);
        }
        return properties;
    }

    /**
     * Forward to Plotly.restyle.
     * Updates specific properties of traces without re-rendering the entire plot
     *
     * @param data properties to update
     * @param traces optional, indices of traces or a single trace index to update
     */
    private _restyle(data: Partial<Data>, traces?: number | number[]) {
        Plotly.restyle(this._plot, data, traces).catch((e: unknown) =>
            setTimeout(() => {
                throw e;
            })
        );
    }

    /**
     * Forward to Plotly.relayout
     * Updates the layout properties of the plot
     *
     * @param layout layout properties to update
     */
    private _relayout(layout: Partial<Layout>) {
        Plotly.relayout(this._plot, layout).catch((e: unknown) =>
            setTimeout(() => {
                throw e;
            })
        );
    }

    /**
     * Forward to Plotly.react
     * Updates the Plotly plot with new data and layout with re-rendering
     *
     * @param traces array of data traces to update
     * @param layout layout properties to update
     */
    private _react(traces: Plotly.Data[], layout: Partial<Layout>): Promise<void> {
        return new Promise<void>((resolve, reject) => {
            Plotly.react(this._plot, traces, layout)
                .then(() => {
                    resolve();
                })
                // eslint-disable-next-line @typescript-eslint/use-unknown-in-catch-callback-variable
                .catch((error: Error) => {
                    setTimeout(() => {
                        reject(error);
                    });
                });
        });
    }

    /** Add all the required callback to the settings */
    private _connectSettings() {
        // Send a warning if a property contains negative values, that will be
        // discarded when using a log scale for this axis
        const negativeLogWarning = (axis: AxisOptions) => {
            if (
                axis.scale.value === 'log' &&
                arrayMaxMin(this._coordinates(axis, 0)[0] as number[])['min'] < 0 &&
                axis.min.value <= 0
            ) {
                this.warnings.sendMessage(
                    'This property contains negative values. Note that taking the log will discard them.'
                );
            }
        };

        // ======= x axis settings
        this._options.x.property.onchange.push(() => {
            negativeLogWarning(this._options.x);
            const values = this._coordinates(this._options.x) as number[][];
            this._restyle({ x: values }, [0, 1]);
            this._relayout({
                'scene.xaxis.title': this._title(this._options.x.property.value),
                'xaxis.title': this._title(this._options.x.property.value),
            } as unknown as Layout);

            if (this._is3D()) {
                this._relayout({
                    'scene.xaxis.autorange': true,
                } as unknown as Layout);
            } else {
                this._relayout({
                    'xaxis.autorange': true,
                } as unknown as Layout);
            }
            this._setScaleStep(this._getBounds().x, 'x');
        });

        this._options.x.scale.onchange.push(() => {
            negativeLogWarning(this._options.x);
            this._options.setLogLabel(this._options.x, 'x');
            if (this._is3D()) {
                this._relayout({
                    'scene.xaxis.type': this._options.x.scale.value,
                } as unknown as Layout);
            } else {
                this._relayout({ 'xaxis.type': this._options.x.scale.value as Plotly.AxisType });
            }
        });

        // function creating a function to be used as onchange callback
        // for <axis>.min and <axis>.max
        const rangeChange = (name: string, axis: AxisOptions, minOrMax: 'min' | 'max') => {
            return (_: number, origin: OptionModificationOrigin) => {
                if (origin === 'JS') {
                    // prevent recursion: this function calls relayout, which then
                    // calls _afterplot, which reset the min/max values.
                    return;
                }
                const min = axis.min.value;
                const max = axis.max.value;
                if (min > max) {
                    this.warnings.sendMessage(
                        `The inserted min and max values in ${name} are such that min > max! The last inserted value was reset.`
                    );
                    if (minOrMax === 'min') {
                        axis.min.reset();
                    } else {
                        axis.max.reset();
                    }
                    return;
                }

                negativeLogWarning(axis);

                if (this._is3D()) {
                    this._relayout({
                        [`scene.${name}.range`]: [min, max],
                    } as unknown as Layout);
                } else {
                    this._relayout({ [`${name}.range`]: [min, max] });
                }
            };
        };

        this._options.x.min.onchange.push(rangeChange('xaxis', this._options.x, 'min'));
        this._options.x.max.onchange.push(rangeChange('xaxis', this._options.x, 'max'));

        // ======= y axis settings
        this._options.y.property.onchange.push(() => {
            negativeLogWarning(this._options.y);
            const values = this._coordinates(this._options.y) as number[][];
            this._restyle({ y: values }, [0, 1]);
            this._relayout({
                'scene.yaxis.title': this._title(this._options.y.property.value),
                'yaxis.title': this._title(this._options.y.property.value),
            } as unknown as Layout);

            if (this._is3D()) {
                this._relayout({
                    'scene.yaxis.autorange': true,
                } as unknown as Layout);
            } else {
                this._relayout({
                    'yaxis.autorange': true,
                } as unknown as Layout);
            }
            this._setScaleStep(this._getBounds().y, 'y');
        });

        this._options.y.scale.onchange.push(() => {
            negativeLogWarning(this._options.y);
            this._options.setLogLabel(this._options.y, 'y');
            if (this._is3D()) {
                this._relayout({
                    'scene.yaxis.type': this._options.y.scale.value,
                } as unknown as Layout);
            } else {
                this._relayout({ 'yaxis.type': this._options.y.scale.value as Plotly.AxisType });
            }
        });

        this._options.y.min.onchange.push(rangeChange('yaxis', this._options.y, 'min'));
        this._options.y.max.onchange.push(rangeChange('yaxis', this._options.y, 'max'));

        // ======= z axis settings
        // setup initial state of the z axis settings
        if (this._options.z.property.value === '') {
            this._options.z.disable();
        } else {
            this._options.z.enable();
        }

        this._options.z.property.onchange.push(() => {
            negativeLogWarning(this._options.z);
            // eslint-disable-next-line @typescript-eslint/no-unsafe-member-access, @typescript-eslint/no-explicit-any
            const was3D = (this._plot as any)._fullData[0].type === 'scatter3d';
            if (this._options.z.property.value === '') {
                if (was3D) {
                    this._switch2D();
                }
                return;
            } else {
                if (!was3D) {
                    this._switch3D();
                }
            }

            const values = this._coordinates(this._options.z);
            this._restyle({ z: values } as Data, [0, 1]);
            this._relayout({
                'scene.zaxis.title': this._title(this._options.z.property.value),
                'scene.zaxis.autorange': true,
            } as unknown as Layout);
            if (this._is3D()) {
                this._setScaleStep(this._getBounds().z as number[], 'z');
            }
        });

        this._options.z.scale.onchange.push(() => {
            negativeLogWarning(this._options.z);
            this._options.setLogLabel(this._options.z, 'z');
            if (this._options.z.property.value !== '') {
                this._relayout({
                    'scene.zaxis.type': this._options.z.scale.value,
                } as unknown as Layout);
            }
        });

        this._options.z.min.onchange.push(rangeChange('zaxis', this._options.z, 'min'));
        this._options.z.max.onchange.push(rangeChange('zaxis', this._options.z, 'max'));

        // ======= color axis settings
        if (this._options.hasColors()) {
            // setup initial color range (must do before setting up events)
            const determineColorRange = (optionMin: number, optionMax: number) => {
                const [min, max] = this._getAxisRange(optionMin, optionMax, 'map.color');
                const minProvided = min !== undefined;
                const maxProvided = max !== undefined;

                const getMinMaxFromValues = () => {
                    const values = this._colors(0)[0] as number[];
                    return arrayMaxMin(values);
                };

                if (minProvided || maxProvided) {
                    // Calculate min/max value from values
                    const { min, max } = getMinMaxFromValues();
                    return {
                        min: minProvided ? optionMin : min,
                        max: maxProvided ? optionMax : max,
                    };
                }

                // Calculate default range
                return getMinMaxFromValues();
            };
            const { min, max } = determineColorRange(
                this._options.color.min.value,
                this._options.color.max.value
            );
            this._options.color.min.value = min;
            this._options.color.max.value = max;
            this._setScaleStep([min, max], 'color');
        } else {
            this._options.color.min.value = 0;
            this._options.color.max.value = 0;
        }

        this._options.color.property.onchange.push(() => {
            if (this._options.hasColors()) {
                this._options.color.mode.enable();
                this._options.color.min.enable();
                this._options.color.max.enable();

                this._colorReset.disabled = false;

                const values = this._colors(0)[0] as number[];
                // Color mode warning needs to be called before setting min and max to avoid isFinite error
                if (canChangeColors(values, 'property')) {
                    const { min, max } = arrayMaxMin(values);
                    // We have to set max first and min second here to avoid sending
                    // a spurious warning in `colorRangeChange` below in case the
                    // new min is bigger than the old max.
                    this._options.color.min.value = Number.NEGATIVE_INFINITY;
                    this._options.color.max.value = max;
                    this._options.color.min.value = min;
                    this._setScaleStep([min, max], 'color');

                    this._relayout({
                        'coloraxis.colorbar.title.text': this._colorTitle(),
                        'coloraxis.showscale': true,
                    } as unknown as Layout);
                }
            } else {
                this._options.color.mode.disable();
                this._options.color.min.disable();
                this._options.color.max.disable();

                this._colorReset.disabled = true;

                this._options.color.min.value = 0;
                this._options.color.max.value = 0;

                this._relayout({
                    'coloraxis.colorbar.title.text': undefined,
                    'coloraxis.showscale': false,
                } as unknown as Layout);
            }

            this._restyle(
                {
                    hovertemplate: this._options.hovertemplate(),
                    'marker.color': this._colors(0),
                    'marker.opacity': this._options.color.opacity.value / 100,
                },
                [0]
            );
        });

        const colorRangeChange = (minOrMax: 'min' | 'max') => {
            const min = this._options.color.min.value;
            const max = this._options.color.max.value;
            if (min > max) {
                this.warnings.sendMessage(
                    `The inserted min and max values in color are such that min > max! The last inserted value was reset.`
                );
                if (minOrMax === 'min') {
                    this._options.color.min.reset();
                } else {
                    this._options.color.max.reset();
                }
                return;
            }

            this._relayout({
                'coloraxis.cmax': max,
                'coloraxis.cmin': min,
                // looks like changing only 'coloraxis.cmax'/'coloraxis.cmin' do
                // not update the color of the points (although it does change
                // the colorbar). Asking for an update of 'coloraxis.colorscale'
                // seems to do the trick. This is possibly a Plotly bug, we
                // would need to investigate a bit more.
                'coloraxis.colorscale': this._options.colorScale(),
            } as unknown as Layout);
        };

        const canChangeColors = (values: number[], changed: string): boolean => {
            const mode = this._options.color.mode.value;

            let invalidValues = '';
            if (mode === 'log' || mode === 'sqrt') {
                invalidValues = '<= 0';
            } else if (mode === 'inverse') {
                invalidValues = '== 0';
            }

            const allValuesNaN = values.every((value) => isNaN(value));
            const someValuesNaN = values.some((value) => isNaN(value));

            if (allValuesNaN) {
                this.warnings.sendMessage(
                    `The selected property contains only values ${invalidValues}. ` +
                        'To display this property, select an appropriate color scale. ' +
                        `The ${changed} will be set to its last value.`
                );

                if (changed === 'property') {
                    this._options.color.property.reset();
                } else {
                    this._options.color.mode.reset();
                }

                return false;
            } else if (someValuesNaN) {
                this.warnings.sendMessage(
                    `The selected property contains some values ${invalidValues}. ` +
                        'These values will be colored in grey.'
                );
                return true;
            } else {
                return true;
            }
        };

        this._options.color.mode.onchange.push(() => {
            const values = this._colors(0)[0] as number[];
            // Color mode warning needs to be called before setting min and max to avoid isFinite error
            if (canChangeColors(values, 'color scale')) {
                const { min, max } = arrayMaxMin(values);
                // We have to set min to infinity first, then max, and then min here
                // to avoid sending a spurious warning in `colorRangeChange` below
                // in case the new min is bigger than the old max.
                this._options.color.min.value = Number.NEGATIVE_INFINITY;
                this._options.color.max.value = max;
                this._options.color.min.value = min;
                this._setScaleStep([min, max], 'color');

                this._relayout({
                    'coloraxis.colorbar.title.text': this._colorTitle(),
                    'coloraxis.showscale': true,
                } as unknown as Layout);

                this._restyle(
                    {
                        hovertemplate: this._options.hovertemplate(),
                        'marker.color': this._colors(0),
                        'marker.opacity': this._options.color.opacity.value / 100,
                    },
                    [0]
                );
            }
        });

        this._options.color.min.onchange.push(() => {
            colorRangeChange('min');
        });
        this._options.color.max.onchange.push(() => {
            colorRangeChange('max');
        });

        this._colorReset.onclick = () => {
            const values = this._colors(0)[0] as number[];
            const { min, max } = arrayMaxMin(values);
            this._options.color.min.value = min;
            this._options.color.max.value = max;
            this._relayout({
                'coloraxis.cmax': max,
                'coloraxis.cmin': min,
                // same as above regarding update of the points color
                'coloraxis.colorscale': this._options.colorScale(),
            } as unknown as Layout);
        };

        // setup initial state of the color GUI
        this._options.color.property.enable();
        if (this._options.hasColors()) {
            this._options.color.mode.enable();
            this._options.color.min.enable();
            this._options.color.max.enable();
            this._options.color.opacity.enable();
            this._colorReset.disabled = false;
        } else {
            this._options.color.min.disable();
            this._options.color.max.disable();
            this._colorReset.disabled = true;
        }

        // ======= color palette
        this._options.color.palette.onchange.push(() => {
            this._relayout({
                'coloraxis.colorscale': this._options.colorScale(),
            } as unknown as Layout);
        });

        // ======= opacity
        this._options.color.opacity.onchange.push(() => {
            this._restyle({
                'marker.opacity': this._options.color.opacity.value / 100,
            });
        });

        // ======= markers symbols
        this._options.symbol.onchange.push(() => {
            this._restyle({ 'marker.symbol': this._symbols() }, [0, 1]);

            this._restyle({
                name: this._legendNames(),
                showlegend: this._showlegend(),
            } as unknown as Data);

            this._relayout({
                'coloraxis.colorbar.len': this._colorbarLen(),
            } as unknown as Layout);
        });

        // ======= markers size
        // setup initial state of the marker size settings
        if (this._options.size.property.value === '') {
            this._options.size.mode.disable();
        } else {
            this._options.size.mode.enable();
        }

        this._options.size.property.onchange.push(() => {
            if (this._options.size.property.value !== '') {
                this._options.size.mode.enable();
            } else {
                this._options.size.mode.disable();
            }
            this._restyle({ 'marker.size': this._sizes(0) } as Data, 0);
        });

        this._options.size.factor.onchange.push(() => {
            this._restyle({ 'marker.size': this._sizes(0) } as Data, 0);
        });

        this._options.size.mode.onchange.push(() => {
            this._restyle({ 'marker.size': this._sizes(0) } as Data, 0);
        });

        this._options.markerOutline.onchange.push(() => {
            const width = this._options.markerOutline.value ? 0.5 : 0;
            this._restyle({ 'marker.line.width': width } as Data, [0]);
        });

        this._options.joinPoints.onchange.push(() => {
            const mode = this._options.joinPoints.value ? 'lines+markers' : 'markers';
            this._restyle({ mode: mode } as Data, [0]);
        });
    }

    /** Actually create the Plotly plot */
    private _createPlot() {
        this._plot.innerHTML = '';

        // Get plot data
        const traces = this._getTraces();

        // Build layout from the options of the settings
        const layout = this._getLayout();

        // Create an empty plot and fill it below
        Plotly.newPlot(this._plot, traces, layout, DEFAULT_CONFIG as unknown as Config)
            .then(() => {
                // In some cases (e.g. in Jupyter notebooks) plotly does not comply
                // with the dimensions of its container unless it receives a resize
                // event _after_ it has loaded. This triggers the event.
                window.requestAnimationFrame(() => {
                    window.dispatchEvent(new Event('resize'));
                });
            })
            .catch((e: unknown) =>
                setTimeout(() => {
                    throw e;
                })
            );
        this._plot.classList.add('chsp-map');
        this._plotFix = fixPlot(this._plot);

        this._plot.on('plotly_click', (event: Plotly.PlotMouseEvent) => {
            // don't update selected env on double click, since it is bound to
            // 'reset zoom level' in 2D mode.
            if (event.event && event.event.detail === 2) {
                return;
            }

            let environment = event.points[0].pointNumber;
            // If clicking on main trace (curveNumber 0) in 2D and LOD is active, remap:
            if (!this._is3D() && event.points[0].curveNumber === 0 && this._lodMap.length) {
                environment = this._lodMap[environment] ?? environment;
            }

            if (this._is3D() && event.points[0].data.name === 'selected') {
                // if someone has clicked on a selection marker, set to active
                // this is only used in 3D mode, since in 2D the HTML marker
                // directly deal with the click event
                for (const [i, [guid, data]] of enumerate(this._selected.entries())) {
                    if (event.points[0].pointNumber === i) {
                        environment = data.current;
                        if (this._active !== guid) {
                            this.setActive(guid);
                            this.activeChanged(
                                guid,
                                this._indexer.fromEnvironment(data.current, this._target)
                            );
                        }
                        break;
                    }
                }
            }

            const indexes = this._indexer.fromEnvironment(environment, this._target);

            this.select(indexes);
            this.onselect(indexes);
        });

        this._plot.on('plotly_relayout', () => this._scheduleLodUpdate());
        this._plot.on('plotly_afterplot', () => this._afterplot());
        this._updateMarkers();

        // set step of min/max select arrows based on the plot range
        const bounds = this._getBounds();
        this._setScaleStep(bounds.x, 'x');
        this._setScaleStep(bounds.y, 'y');
        if (bounds.z !== undefined) {
            this._setScaleStep(bounds.z, 'z');
        }

        // Hack to fix a Plotly bug preventing zooming on Safari
        this._plot.addEventListener('wheel', () => {});
    }

    /**
     * Builds the layout to be provided to Plotly from the options
     */
    private _getLayout(): Partial<Layout> {
        // make a copy of the default layout
        const layout = JSON.parse(JSON.stringify(DEFAULT_LAYOUT)) as typeof DEFAULT_LAYOUT;
        // and set values specific to the displayed dataset
        layout.xaxis.title = this._title(this._options.x.property.value);
        layout.yaxis.title = this._title(this._options.y.property.value);
        layout.xaxis.type = this._options.x.scale.value;
        layout.yaxis.type = this._options.y.scale.value;
        layout.scene.xaxis.title = this._title(this._options.x.property.value);
        layout.scene.yaxis.title = this._title(this._options.y.property.value);
        layout.scene.zaxis.title = this._title(this._options.z.property.value);
        layout.coloraxis.colorscale = this._options.colorScale();
        layout.coloraxis.cmin = this._options.color.min.value;
        layout.coloraxis.cmax = this._options.color.max.value;
        layout.coloraxis.colorbar.title.text = this._colorTitle();
        layout.coloraxis.colorbar.len = this._colorbarLen();
        layout.coloraxis.showscale = this._options.hasColors();

        // Set ranges for the axes
        layout.xaxis.range = this._getAxisRange(
            this._options.x.min.value,
            this._options.x.max.value,
            'map.x'
        );
        layout.yaxis.range = this._getAxisRange(
            this._options.y.min.value,
            this._options.y.max.value,
            'map.y'
        );
        layout.zaxis.range = this._getAxisRange(
            this._options.z.min.value,
            this._options.z.max.value,
            'map.z'
        );
        return layout as Partial<Layout>;
    }

    /**
     * Validate min/max options provided by the user. We use `NaN` internally to mark missing values,
     * which are then transformed into undefined by this function
     */
    private _getAxisRange = (
        min: number,
        max: number,
        axisName: string
    ): [number | undefined, number | undefined] => {
        const minProvided = !isNaN(min);
        const maxProvided = !isNaN(max);

        // At least one range value is specified. By default, zeros are set
        if (minProvided && maxProvided) {
            if (min <= max) {
                return [min, max];
            }
            this.warnings.sendMessage(
                `The inserted min and max values in ${axisName} are such that min > max!` +
                    `The default values will be used.`
            );
        }
        return [minProvided ? min : undefined, maxProvided ? max : undefined];
    };

    /** Get the property with the given name */
    private _property(name: string): NumericProperty {
        const result = this._data[this._target][name];
        if (result === undefined) {
            throw Error(`unknown property '${name}' requested in map`);
        }
        return result;
    }

    /**
     * Get the values associated with the given `axis`, to use with the given
     * plotly `trace`, or all of them if `trace === undefined`
     *
     * @param  axis   Options of the axis we need coordinates for
     * @param  trace  plotly trace for which we require coordinate
     * @return        data usable with Plotly.restyle
     */
    private _coordinates(axis: AxisOptions, trace?: number): Array<undefined | number[]> {
        // this happen for the z axis in 2D mode
        if (axis.property.value === '') {
            return this._selectTrace(undefined, undefined, trace);
        }

        const values = this._property(axis.property.value).values;
        // in 2d mode, set all selected markers coordinates to NaN since we are
        // using HTML markers instead.
        const selected = [];
        for (const marker of this._selected.values()) {
            if (this._is3D()) {
                selected.push(values[marker.current]);
            } else {
                selected.push(NaN);
            }
        }
        return this._selectTrace<number[]>(values, selected, trace);
    }

    private _title(name: string): string {
        let propertyTitle = name;
        if (name !== '') {
            const units = this._property(name).units;
            if (units !== undefined) {
                propertyTitle = name + ` / ${units}`;
            }
        }
        return propertyTitle;
    }

    private _colorTitle(): string {
        let title = this._title(this._options.color.property.value);
        switch (this._options.color.mode.value) {
            case 'inverse':
                title = `(${title})<sup>-1</sup>`;
                break;
            case 'log':
                title = `log<sub>10</sub>(${title})`;
                break;
            case 'sqrt':
                title = `&#x221A;(${title})`;
                break;
            case 'linear':
                break;
            default:
                break;
        }
        return title;
    }

    /**
     * Get the color values to use with the given plotly `trace`, or all of
     * them if `trace === undefined`
     */
    private _colors(trace?: number): Array<Array<string | number>> {
        let colors;
        if (this._options.hasColors()) {
            colors = this._property(this._options.color.property.value).values;
        } else {
            colors = new Array(this._property(this._options.x.property.value).values.length).fill(
                0.5
            ) as number[];
        }
        const values = this._options.calculateColors(colors);
        const selected = [];
        for (const data of this._selected.values()) {
            selected.push(data.color);
        }

        return this._selectTrace<Array<string | number>>(values, selected, trace);
    }

    /**
     * Get the values to use as marker size with the given plotly `trace`, or
     * all of them if `trace === undefined`.
     */
    private _sizes(trace?: number): Array<number | number[]> {
        let sizes;
        if (this._options.size.property.value !== '') {
            sizes = this._property(this._options.size.property.value).values;
        } else {
            sizes = new Array(this._property(this._options.x.property.value).values.length).fill(
                1.0
            ) as number[];
        }
        const values = this._options.calculateSizes(sizes);
        const selected = [];
        if (this._is3D()) {
            for (const guid of this._selected.keys()) {
                if (guid === this._active) {
                    selected.push(1000);
                } else {
                    selected.push(500);
                }
            }
        }
        return this._selectTrace<number | number[]>(values, selected, trace);
    }

    /**
     * Get the values to use as marker symbol with the given plotly `trace`, or
     * all of them if `trace === undefined`.
     */
    private _symbols(trace?: number): Array<string | string[] | number[]> {
        if (this._options.symbol.value === '') {
            // default to 0 (i.e. circles)
            return this._selectTrace<string | string[]>('circle', 'circle', trace);
        }

        const property = this._property(this._options.symbol.value);
        const symbols = this._options.getSymbols(property);
        const selected = [];
        for (const data of this._selected.values()) {
            selected.push(symbols[data.current]);
        }
        return this._selectTrace<typeof symbols>(symbols, selected as typeof symbols, trace);
    }

    /** Should we show the legend for the various symbols used? */
    private _showlegend(): boolean[] {
        const result = [false, false];

        if (this._options.symbol.value !== '') {
            for (let i = 0; i < this._symbolsCount(); i++) {
                result.push(true);
            }
            return result;
        } else {
            for (let i = 0; i < this._data.maxSymbols; i++) {
                result.push(false);
            }
            return result;
        }
    }

    /** Get the list of symbols names to use for the legend */
    private _legendNames(): string[] {
        const result = ['', ''];

        if (this._options.symbol.value !== '') {
            const property = this._property(this._options.symbol.value);
            assert(property.string !== undefined);
            for (const name of property.string.strings()) {
                result.push(name);
            }
            return result;
        } else {
            for (let i = 0; i < this._data.maxSymbols; i++) {
                result.push('');
            }
            return result;
        }
    }

    /**
     * Select either main, selected or both depending on `trace`, and return
     * them in a mode usable with `Plotly.restyle`/{@link PropertiesMap._restyle}
     */
    private _selectTrace<T>(main: T, selected: T, trace?: number): T[] {
        if (trace === 0) {
            return [main];
        } else if (trace === 1) {
            return [selected];
        } else if (trace === undefined) {
            return [main, selected];
        } else {
            throw Error('internal error: invalid trace number');
        }
    }

    /** Get the length of the colorbar to accommodate for the legend */
    private _colorbarLen(): number {
        /// Heigh of a legend item in plot unit
        const LEGEND_ITEM_HEIGH = 0.045;
        const PADDING = 0.025;
        return 1 - LEGEND_ITEM_HEIGH * this._symbolsCount() - PADDING;
    }

    /** Is the the current plot a 3D plot? */
    private _is3D(): boolean {
        return this._options.is3D();
    }

    /** How many symbols are on this plot?*/
    private _symbolsCount(): number {
        if (this._options.symbol.value !== '') {
            const property = this._property(this._options.symbol.value);
            assert(property.string !== undefined);
            return property.string.strings().length;
        } else {
            return 0;
        }
    }

    /** Switch current plot from 2D to 3D */
    private _switch3D(): void {
        assert(this._is3D());
        this._options.z.enable();

        const symbols = this._symbols();
        for (let s = 0; s < this._data.maxSymbols; s++) {
            symbols.push([get3DSymbol(s)]);
        }

        // switch all traces to 3D mode
        this._restyle({
            'marker.symbol': symbols,
            type: 'scatter3d',
        } as unknown as Data);

        for (const data of this._selected.values()) {
            data.toggleVisible(false);
        }
        this._updateMarkers();

        // Change the data that vary between 2D and 3D mode
        const marker_line = this._options.markerOutline.value ? 0.5 : 0.0;
        this._restyle(
            {
                'marker.line.width': marker_line,
                // size change from 2D to 3D
                'marker.size': this._sizes(),
                'marker.sizemode': 'area',
            } as Data,
            [0, 1]
        );

        this._relayout({
            // change colorbar length to accommodate for symbols legend
            'coloraxis.colorbar.len': this._colorbarLen(),
            // Carry over axis types
            'scene.xaxis.type': this._options.x.scale.value as Plotly.AxisType,
            'scene.yaxis.type': this._options.y.scale.value as Plotly.AxisType,
            'scene.zaxis.type': this._options.z.scale.value as Plotly.AxisType,
        } as unknown as Layout);
    }

    /** Switch current plot from 3D back to 2D */
    private _switch2D(): void {
        assert(!this._is3D());
        this._options.z.disable();

        const symbols = this._symbols();
        for (let sym = 0; sym < this._data.maxSymbols; sym++) {
            symbols.push([sym]);
        }

        // switch all traces to 2D mode
        this._restyle({
            'marker.symbol': symbols,
            type: 'scattergl',
        } as unknown as Data);

        // show selected environments markers
        for (const data of this._selected.values()) {
            data.toggleVisible(true);
        }

        const marker_line = this._options.markerOutline.value ? 0.5 : 0.0;
        this._restyle(
            {
                'marker.line.width': marker_line,
                // size change from 3D to 2D
                'marker.size': this._sizes(),
            } as Data,
            [0, 1]
        );

        this._relayout({
            // change colorbar length to accommodate for symbols legend
            'coloraxis.colorbar.len': this._colorbarLen(),
            // Carry over axis types
            'xaxis.type': this._options.x.scale.value as Plotly.AxisType,
            'yaxis.type': this._options.y.scale.value as Plotly.AxisType,
        } as unknown as Layout);
    }

    /**
     * Function used as callback to update the axis ranges in settings after
     * the user changes zoom or range on the plot
     */
    private _afterplot(): void {
        const bounds = this._getBounds();
        const updateAxisValues = (axis: AxisOptions, [boundMin, boundMax]: [number, number]) => {
            axis.min.value = isNaN(axis.min.value) ? boundMin : axis.min.value;
            axis.max.value = isNaN(axis.max.value) ? boundMax : axis.max.value;
        };

        // Set calculated value to the range if axis bound value was not provided by user
        updateAxisValues(this._options.x, bounds.x);
        updateAxisValues(this._options.y, bounds.y);
        if (bounds.z !== undefined) {
            updateAxisValues(this._options.z, bounds.z);
        }

        if (!this._is3D()) {
            this._updateMarkers();
        }
    }

    /**
     * Update the position, color & size of markers. If `markers` is present, only
     * markers inside the array are updated, otherwise everything is updated.
     */
    private _updateMarkers(markers?: MarkerData[]): void {
        if (markers === undefined) {
            markers = Array.from(this._selected.values());
        }

        if (this._is3D()) {
            markers.forEach((marker) => marker.toggleVisible(false));
            this._restyle(
                {
                    'marker.color': this._colors(1),
                    'marker.opacity': this._options.color.opacity.value / 100,
                    'marker.size': this._sizes(1),
                    'marker.symbol': this._symbols(1),
                    x: this._coordinates(this._options.x, 1),
                    y: this._coordinates(this._options.y, 1),
                    z: this._coordinates(this._options.z, 1),
                } as Data,
                1
            );
        } else {
            const allX = this._coordinates(this._options.x, 0)[0] as number[];
            const allY = this._coordinates(this._options.y, 0)[0] as number[];
            const plotWidth = this._plot.getBoundingClientRect().width;

            for (const marker of markers) {
                let x = allX[marker.current];
                let y = allY[marker.current];

                // make sure we deal with log scale
                if (this._options.x.scale.value === 'log') {
                    x = Math.log10(x);
                }
                if (this._options.y.scale.value === 'log') {
                    y = Math.log10(y);
                }

                // convert to pixel coordinates
                x = plotWidth - this._pixelCoordinate(x, 'x');
                y = this._pixelCoordinate(y, 'y');

                const bounds = this._getBounds();
                const xMin = plotWidth - this._pixelCoordinate(bounds.x[1], 'x');
                const xMax = plotWidth - this._pixelCoordinate(bounds.x[0], 'x');

                // pixel coordinates are inverted for y (y is measured from the
                // top of the page, not the bottom)
                const yMin = this._pixelCoordinate(bounds.y[1], 'y');
                const yMax = this._pixelCoordinate(bounds.y[0], 'y');

                const isInsideRange = (value: number, min: number, max: number) => {
                    // allow points a bit outside of the range, according to
                    // this the tolerance value in pixels
                    const tolerance = 10;
                    return value + tolerance > min && value - tolerance < max;
                };

                if (isInsideRange(x, xMin, xMax) && isInsideRange(y, yMin, yMax)) {
                    marker.toggleVisible(true);
                    marker.update({ x, y });
                } else {
                    marker.toggleVisible(false);
                }
            }
        }
    }

    // Get the current boundaries on x/y/z axis
    private _getBounds(): { x: [number, number]; y: [number, number]; z?: [number, number] } {
        if (this._is3D()) {
            // HACK: `_fullLayout` is not public, so it might break
            const layout = this._plot._fullLayout.scene;
            return {
                x: layout.xaxis.range as [number, number],
                y: layout.yaxis.range as [number, number],
                z: layout.zaxis.range as [number, number],
            };
        } else {
            const layout = this._plot._fullLayout;
            return {
                x: layout.xaxis.range as [number, number],
                y: layout.yaxis.range as [number, number],
            };
        }
    }

    // Computes the pixel coordinate inside the plot div of a given value along
    // the given axis
    private _pixelCoordinate(value: number, axisName: 'x' | 'y'): number {
        assert(!this._is3D());
        let axis;
        switch (axisName) {
            case 'x':
                axis = this._plot._fullLayout.xaxis;
                break;
            case 'y':
                axis = this._plot._fullLayout.yaxis;
                break;
        }

        return axis.l2p(value) + axis._offset;
    }

    /** Changes the step of the arrow buttons in min/max input based on dataset range*/
    private _setScaleStep(axisBounds: number[], name: 'x' | 'y' | 'z' | 'color'): void {
        if (axisBounds !== undefined) {
            // round to 10 decimal places so it does not break in Firefox
            const step = Math.round(((axisBounds[1] - axisBounds[0]) / 20) * 10 ** 10) / 10 ** 10;
            const minElement = this._options.getModalElement<HTMLInputElement>(`map-${name}-min`);
            const maxElement = this._options.getModalElement<HTMLInputElement>(`map-${name}-max`);
            minElement.step = `${step}`;
            maxElement.step = `${step}`;
        }
    }
}

/** Extract the data associated with the first `path` element in an SVG string */
function extractSvgPath(svg: string) {
    const doc = document.createElement('div');
    doc.innerHTML = svg;
    return doc.getElementsByTagName('path')[0].getAttribute('d');
}
