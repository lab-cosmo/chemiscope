/**
 * @packageDocumentation
 * @module map
 */

import assert from 'assert';

import Plotly from './plotly/plotly-scatter';
import { Config, Data, Layout, PlotlyScatterElement } from './plotly/plotly-scatter';
import * as plotlyStyles from './plotly/plotly-styles';

import { Property, Settings } from '../dataset';

import { DisplayTarget, EnvironmentIndexer, Indexes } from '../indexer';
import { OptionModificationOrigin } from '../options';
import { cameraToPlotly, plotlyToCamera } from '../utils/camera';
import { Bounds, GUID, PositioningCallback, Warnings, arrayMaxMin } from '../utils';
import { enumerate, getElement, getFirstKey } from '../utils';

import { MapData, NumericProperties, NumericProperty } from './data';
import { MarkerData } from './marker';
import { AxisOptions, MapOptions, get3DSymbol } from './options';
import { computeLODIndices, computeScreenSpaceLOD } from './lod';
import * as styles from '../styles';

import PNG_SVG from '../static/download-png.svg';
import SVG_SVG from '../static/download-svg.svg';

/**
 * Export the plot as a PNG or SVG image, hiding the "selected" trace (index 1)
 * during the export.
 *
 * @param gd Plotly plot element
 * @param format format of the image
 */
function exportImage(gd: PlotlyScatterElement, format: 'png' | 'svg') {
    const width = Math.max(gd._fullLayout.width, 600);
    const ratio = gd._fullLayout.height / gd._fullLayout.width;

    const opts: Plotly.DownloadImgopts = {
        filename: 'chemiscope-map',
        format: format,
        width: width,
        height: format === 'png' ? width * ratio : Math.max(gd._fullLayout.height, 600),
    };

    if (format === 'png') {
        // scale is not part of `DownloadImgopts`, but accepted
        // by the function anyway
        (opts as unknown as { scale: number }).scale = 3;
    }

    // Hide the "selected" trace (index 1) for the export.
    // In 2D mode, this trace is already empty (using NaNs), but in 3D mode
    // it contains the markers for selected environments.
    Plotly.restyle(gd, { visible: false }, [1])
        .then(() => {
            return Plotly.downloadImage(gd, opts);
        })
        .then(() => {
            return Plotly.restyle(gd, { visible: true }, [1]);
        })
        .catch((e: unknown) => {
            // make sure we show the trace again even if download failed
            void Plotly.restyle(gd, { visible: true }, [1]);
            setTimeout(() => {
                throw e;
            });
        });
}

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
    dragmode: 'zoom',
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
            autorange: true,
            range: undefined as (number | undefined)[] | undefined,
            title: { text: '' },
            type: 'linear',
        },
        yaxis: {
            showspikes: false,
            autorange: true,
            range: undefined as (number | undefined)[] | undefined,
            title: { text: '' },
            type: 'linear',
        },
        zaxis: {
            showspikes: false,
            autorange: true,
            range: undefined as (number | undefined)[] | undefined,
            title: { text: '' as undefined | string },
            type: 'linear',
        },
    },
    showlegend: true,
    xaxis: {
        range: undefined as (number | undefined)[] | undefined,
        title: { text: '' },
        type: 'linear',
        zeroline: false,
    },
    yaxis: {
        range: undefined as (number | undefined)[] | undefined,
        title: { text: '' },
        type: 'linear',
        zeroline: false,
    },
    zaxis: {
        range: undefined as (number | undefined)[] | undefined,
        title: { text: '' },
        type: 'linear',
        zeroline: false,
    },
};

const DEFAULT_CONFIG = {
    displayModeBar: true,
    displaylogo: false,
    responsive: true,
    doubleClick: 'reset',
    doubleClickDelay: 600,
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
        'zoom3d',
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
                    exportImage(gd, 'png');
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
                    exportImage(gd, 'svg');
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

    /**
     * LOD (Level of Detail) Configuration
     *
     * Speeds up rendering of large datasets by downsampling points
     * when zoomed out.
     */
    private static readonly LOD_THRESHOLD = 50000;
    /// Stores the subset of point indices to display when LOD is active
    private _lodIndices: number[] | null = null;
    /**
     * Guard to skip concurrent LOD updates. When set, other callers simply return early
     * instead of waiting
     */
    private _lodBusy = false;
    // Timeout id used to batch plotly afterplot events
    private _afterplotRequest: number | null = null;

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

        // Determine whether to show the LOD option based on dataset size
        const nPoints = Object.values(currentProperties)[0].values.length;
        if (nPoints > PropertiesMap.LOD_THRESHOLD) {
            this._options.showLODOption(true);
        }

        // Initial LOD computation
        this._computeLOD();

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
    public switchTarget(target: DisplayTarget): Promise<void> {
        // Check if the target value actually changed
        if (target !== this._target) {
            // Set new widget target
            this._target = target;
            // Reset LOD
            this._lodIndices = null;

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

            // Recompute LOD with new options
            this._computeLOD();

            // Append the callbacks to the new options
            this._connectSettings();

            // Re-render the plot with the new data and layout
            return this._react(this._getTraces(), this._getLayout(), this._getConfig());
        }
        return Promise.resolve();
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
     * Export the plot as a PNG data URI, hiding the "selected" trace (index 1)
     * during the export.
     */
    public async exportPNG(): Promise<string> {
        const fullLayout = this._plot._fullLayout;
        const width = Math.max(fullLayout.width, 600);
        const ratio = fullLayout.height / fullLayout.width;
        const height = width * ratio;

        const opts: Plotly.ToImgopts = {
            format: 'png',
            width: width,
            height: height,
            scale: 3,
        };

        // Hide the "selected" trace (index 1) for the export.
        try {
            await Plotly.restyle(this._plot, { visible: false }, [1]);
            const dataUrl = await Plotly.toImage(this._plot, opts);
            await Plotly.restyle(this._plot, { visible: true }, [1]);
            return dataUrl;
        } catch (e) {
            void Plotly.restyle(this._plot, { visible: true }, [1]);
            throw e;
        }
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
            customdata: this._colorValues(0)[0],
            marker: {
                color: this._colors(0)[0],
                coloraxis: 'coloraxis',
                showscale: false,
                line: {
                    color: this._lineColors(0)[0],
                    width: this._options.markerOutline.value ? 0.5 : 0,
                },
                opacity: this._options.color.opacity.value / 100,
                size: this._sizes(0)[0],
                sizemode: 'area',
                symbol: this._symbols(0)[0],
            },
            line: {
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

        // Dummy trace to display the colorbar regardless of the styling of
        // the main trace. Useful when activating selection mode
        const range = this._getColorRange();
        const dummy = {
            name: 'colorbar-dummy',
            type: type,
            x: this._coordinates(this._options.x, 2)[0],
            y: this._coordinates(this._options.y, 2)[0],
            z: this._coordinates(this._options.z, 2)[0],
            marker: {
                color: [range.min, range.max],
                coloraxis: 'coloraxis',
                showscale: true,
                opacity: 0,
                size: 0,
            },
            visible: this._options.hasColors(),
            showlegend: false,
            hoverinfo: 'none',
        };

        const traces = [main as Data, selected as Data, dummy as Data];

        // Calculate legend names and show legend flags based on data properties
        const legendNames = this._legendNames().slice(3);
        const showlegend = this._showlegend().slice(3);
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

                // Update reference to the color reset button, as the old one was removed from DOM
                this._colorReset =
                    this._options.getModalElement<HTMLButtonElement>('map-color-reset');

                // Update LOD toggle visibility based on the new target's dataset size.
                const nPoints = Object.values(properties)[0].values.length;
                if (nPoints > PropertiesMap.LOD_THRESHOLD) {
                    this._options.showLODOption(true);
                }
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
     * @param config optional plot configuration
     */
    private _react(
        traces: Plotly.Data[],
        layout: Partial<Layout>,
        config?: Partial<Config>
    ): Promise<void> {
        return (
            Plotly.react(this._plot, traces, layout, config as Config)
                .then(() => {})
                // eslint-disable-next-line @typescript-eslint/use-unknown-in-catch-callback-variable
                .catch((error: Error) => {
                    setTimeout(() => {
                        throw error;
                    });
                })
        );
    }

    /** Add all the required callback to the settings */
    private _connectSettings() {
        // Range reset button
        const resetRanges = this._options.getModalElement<HTMLButtonElement>('map-range-reset');
        resetRanges.onclick = () => {
            this._resetToGlobalView();
        };

        // Send a warning if a property contains negative values, that will be
        // discarded when using a log scale for this axis
        const negativeLogWarning = (axis: AxisOptions) => {
            // skip if axis has no property set (e.g. z axis in 2d mode)
            if (axis.property.value === '') {
                return;
            }

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
            // Reset min/max when changing property to trigger autoscale
            negativeLogWarning(this._options.x);
            this._handleAxisPropertyChange('x');
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
                        [`scene.${name}.autorange`]: false,
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
            this._handleAxisPropertyChange('y');
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
            // Reset min/max when changing property to trigger autoscale
            this._options.z.min.value = NaN;
            this._options.z.max.value = NaN;
            negativeLogWarning(this._options.z);

            const was3D = this._plot._fullData[0].type === 'scatter3d';

            // If no z property selected -> switch to 2D
            if (this._options.z.property.value === '') {
                if (was3D) {
                    this._switch2D();
                    // ... and re-update LOD in 2D
                    this._updateLOD(this._getBounds());
                }
                return;
            } else {
                if (!was3D) {
                    this._switch3D();
                }
            }

            // LOD: Z changed, compute a first downsampling if necessary
            this._computeLOD();

            Plotly.relayout(this._plot, {
                'scene.zaxis.title.text': this._title(this._options.z.property.value),
                'scene.zaxis.autorange': true,
            } as unknown as Layout)
                .then(() => {
                    // The zrange is now known, and we can trigger a proper subsampling
                    const zRange = this._plot._fullLayout.scene.zaxis.range as number[];
                    this._options.z.min.value = zRange[0];
                    this._options.z.max.value = zRange[1];

                    if (this._is3D()) {
                        this._setScaleStep(this._getBounds().z as number[], 'z');
                    }

                    // re-update LOD based on known ranges
                    this._updateLOD(this._getBounds());
                })
                .catch((e: unknown) => {
                    setTimeout(() => {
                        throw e;
                    });
                });
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
                    return arrayMaxMin(this._getNumericColors());
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

                const numericValues = this._getNumericColors();
                // Color mode warning needs to be called before setting min and max to avoid isFinite error
                if (canChangeColors(numericValues, 'property')) {
                    const { min, max } = arrayMaxMin(numericValues);
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

            this._restyleFull();
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
            const numericValues = this._getNumericColors();
            // Color mode warning needs to be called before setting min and max to avoid isFinite error
            if (canChangeColors(numericValues, 'color scale')) {
                const { min, max } = arrayMaxMin(numericValues);
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

                this._restyleFull();
            }
        });

        this._options.color.min.onchange.push(() => {
            colorRangeChange('min');
        });
        this._options.color.max.onchange.push(() => {
            colorRangeChange('max');
        });

        this._colorReset.onclick = () => {
            const numericValues = this._getNumericColors();
            const { min, max } = arrayMaxMin(numericValues);
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
            this._restyleFull();
        });

        // ======= opacity
        this._options.color.opacity.onchange.push(() => {
            this._restyle({
                'marker.opacity': this._options.color.opacity.value / 100,
            });
        });

        // ======= selection
        const updateColors = () => {
            this._restyleFull();
        };

        this._options.color.select.mode.onchange.push(updateColors);
        this._options.color.select.category.onchange.push(updateColors);
        this._options.color.select.min.onchange.push(updateColors);
        this._options.color.select.max.onchange.push(updateColors);

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
            this._restyle({ 'marker.line.width': width } as Data, 0);
        });

        this._options.joinPoints.onchange.push(() => {
            const mode = this._options.joinPoints.value ? 'lines+markers' : 'markers';
            this._restyle({ mode: mode } as Data, [0]);
        });

        // ======= Level of Detail settings
        this._options.useLOD.onchange.push(() => {
            const bounds = this._getBounds();
            this._computeLOD(bounds);
            // Force a full restyle. Since _lodIndices will be null if disabled,
            // this will render all points.
            this._restyleFull();
        });

        // ======= camera state update
        this._options.camera.onchange.push((camera, origin) => {
            if (origin === 'JS' && camera !== undefined && this._is3D()) {
                const update = cameraToPlotly(camera);
                this._relayout({
                    'scene.camera': update.camera,
                    'scene.aspectratio': update.aspectratio,
                } as unknown as Layout);
            }
        });
    }

    /** Actually create the Plotly plot */
    private _createPlot() {
        this._plot.innerHTML = '';

        // Get plot data
        const traces = this._getTraces();

        // Build layout from the options of the settings
        const layout = this._getLayout();

        // Get config
        const config = this._getConfig();

        // Create an empty plot and fill it below
        Plotly.newPlot(this._plot, traces, layout, config)
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

        // Add callbacks

        this._plot.on('plotly_click', (eventData: Plotly.PlotMouseEvent) => {
            const event = eventData;
            // don't update selected env on double click, since it is bound to
            // 'reset zoom level' in 2D mode.
            if (event.event && event.event.detail === 2) {
                return;
            }

            // don't intercept clicks while the subsampling is updating
            if (this._lodBusy) {
                return;
            }
            let environment = event.points[0].pointNumber;

            // LOD: Map the clicked point back to real environment index if LOD is active on main trace
            if (this._lodIndices !== null && event.points[0].curveNumber === 0) {
                if (environment >= this._lodIndices.length) {
                    // ignore stray clicks that happen during redrawing
                    return;
                }
                environment = this._lodIndices[environment];
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

        this._plot.on('plotly_afterplot', () => {
            if (this._lodBusy) {
                return;
            }
            if (this._afterplotRequest !== null) {
                window.clearTimeout(this._afterplotRequest);
            }

            this._afterplotRequest = window.setTimeout(() => {
                this._afterplotRequest = null;
                this._afterplot();
            }, 50);
        });

        // 3D LOD: Listen to relayout to catch 3D camera changes (zoom/pan)
        this._plot.on('plotly_relayout', (event: Plotly.PlotRelayoutEvent) => {
            if (this._lodBusy) {
                return;
            }
            if (this._afterplotRequest !== null) {
                window.clearTimeout(this._afterplotRequest);
            }
            this._afterplotRequest = window.setTimeout(() => {
                this._afterplotRequest = null;
                this._afterplot();

                // Check if LOD update is needed based on relayout event
                const lodEnabled =
                    this._options.useLOD.value &&
                    this._property(this._options.x.property.value).values.length >
                        PropertiesMap.LOD_THRESHOLD;
                const viewChanged = Object.keys(event).some(
                    (key) =>
                        key.match(/^(xaxis|yaxis)\.range/) ||
                        key.match(/^scene\.(camera|aspectratio)/) ||
                        key.includes('autorange')
                );
                if (lodEnabled && viewChanged) {
                    this._updateLOD(this._getBounds());
                }
            }, 50);
        });

        // Handle double-click to reset view (global LOD)
        this._plot.on('plotly_doubleclick', () => {
            this._resetToGlobalView();
            return false;
        });

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

        // Hack to ensure that the active trace is selected AFTER the plot is drawn
        // which seems to be necessary to ensure _fullLayout actually contains
        // the state of the plotly viewer in 3D
        setTimeout(() => {
            if (this._active !== undefined) {
                this.setActive(this._active);
                const data = this._selected.get(this._active);
                if (data !== undefined) {
                    this.activeChanged(
                        this._active,
                        this._indexer.fromEnvironment(data.current, this._target)
                    );
                }
            }
        }, 1000);
    }

    private _handleAxisPropertyChange(axis: 'x' | 'y'): void {
        const axisOptions = this._options[axis];

        // Reset min/max when changing property to trigger autoscale
        axisOptions.min.value = NaN;
        axisOptions.max.value = NaN;

        // Recompute LOD since axis property changed
        this._computeLOD();

        // Request a full restyle of traces
        this._restyleFull();

        this._relayout({
            [`scene.${axis}axis.title.text`]: this._title(axisOptions.property.value),
            [`${axis}axis.title.text`]: this._title(axisOptions.property.value),
        } as unknown as Layout);

        if (this._is3D()) {
            this._relayout({ [`scene.${axis}axis.autorange`]: true } as unknown as Layout);
        } else {
            this._relayout({ [`${axis}axis.autorange`]: true } as unknown as Layout);
        }

        const bounds = this._getBounds();
        this._setScaleStep(bounds[axis], axis);
    }

    /**
     * Builds the layout to be provided to Plotly from the options
     */
    private _getLayout(): Partial<Layout> {
        // make a copy of the default layout
        const layout = JSON.parse(JSON.stringify(DEFAULT_LAYOUT)) as typeof DEFAULT_LAYOUT;
        // and set values specific to the displayed dataset
        layout.xaxis.title.text = this._title(this._options.x.property.value);
        layout.yaxis.title.text = this._title(this._options.y.property.value);
        layout.xaxis.type = this._options.x.scale.value;
        layout.yaxis.type = this._options.y.scale.value;
        layout.scene.xaxis.title.text = this._title(this._options.x.property.value);
        layout.scene.yaxis.title.text = this._title(this._options.y.property.value);
        layout.scene.zaxis.title.text = this._title(this._options.z.property.value);
        layout.scene.xaxis.type = this._options.x.scale.value;
        layout.scene.yaxis.type = this._options.y.scale.value;
        layout.scene.zaxis.type = this._options.z.scale.value;
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

        // 3D scenes have a separate axis configuration
        if (this._is3D()) {
            layout.scene.xaxis.range = this._getAxisRange(
                this._options.x.min.value,
                this._options.x.max.value,
                'map.x'
            );
            layout.scene.xaxis.autorange = this._getAxisAutoRange(
                this._options.x.min.value,
                this._options.x.max.value
            );
            layout.scene.yaxis.range = this._getAxisRange(
                this._options.y.min.value,
                this._options.y.max.value,
                'map.y'
            );
            layout.scene.yaxis.autorange = this._getAxisAutoRange(
                this._options.y.min.value,
                this._options.y.max.value
            );
            layout.scene.zaxis.range = this._getAxisRange(
                this._options.z.min.value,
                this._options.z.max.value,
                'map.z'
            );
            layout.scene.zaxis.autorange = this._getAxisAutoRange(
                this._options.z.min.value,
                this._options.z.max.value
            );
            layout.dragmode = 'orbit';

            if (this._options.camera.value) {
                const update = cameraToPlotly(this._options.camera.value);
                const scene = layout.scene;
                Object.assign(scene.camera, update.camera);
                Object.assign(scene, { aspectratio: update.aspectratio });
            }
        }

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

    private _getAxisAutoRange = (min: number, max: number): boolean => {
        const minProvided = !isNaN(min);
        const maxProvided = !isNaN(max);

        // Both ranges are provided
        return !(minProvided && maxProvided);
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
            return this._selectTrace(undefined, undefined, undefined, trace);
        }

        const values = this._property(axis.property.value).values;

        // LOD: Apply binning filter for main trace (trace 0)
        // If trace is undefined (both), _selectTrace distributes the first argument to
        // main trace
        let mainValues = values;
        if (this._options.color.select.mode.value.endsWith('hide')) {
            const mask = this._getSelectionMask();
            mainValues = values.map((v, i) => (mask[i] ? v : NaN));
        }

        if (trace === 0 || trace === undefined) {
            mainValues = this._applyLOD(mainValues);
        }

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

        const val0 = values[0] || 0;
        const dummy = [val0, val0];

        return this._selectTrace(mainValues, selected, dummy, trace);
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
     * Get the numeric color values after applying all the scaling and transforms
     * selected by the options.
     */
    private _getNumericColors(): number[] {
        const propValues = this._property(this._options.color.property.value).values;
        const values = this._options.calculateColors(propValues);
        return values.filter((v): v is number => typeof v === 'number');
    }

    private _getColorRange(): { min: number; max: number } {
        let min = this._options.color.min.value;
        let max = this._options.color.max.value;

        if (isNaN(min) || isNaN(max)) {
            const numericValues = this._getNumericColors();
            const { min: dMin, max: dMax } = arrayMaxMin(numericValues);
            if (isNaN(min)) min = dMin;
            if (isNaN(max)) max = dMax;
        }
        return { min, max };
    }

    /**
     * Get the mask of selected points based on the current selection mode
     */
    private _getSelectionMask(): boolean[] {
        const n = this._property(this._options.x.property.value).values.length;
        if (this._options.color.select.mode.value === 'all') {
            return new Array<boolean>(n).fill(true);
        }

        let colors;
        if (this._options.hasColors()) {
            colors = this._property(this._options.color.property.value).values;
        } else {
            // "Fixed" color corresponds to a constant value
            // (should not be used for filtering but should never happen)
            colors = new Array<number>(n).fill(0.5);
        }

        let categoryValues: string[] | undefined;
        const catProp = this._options.getCategorySelectionProperty();
        if (catProp) {
            const prop = this._property(catProp);
            if (prop.string) {
                const interner = prop.string;
                categoryValues = prop.values.map((v) => interner.string(v));
            }
        }
        return this._options.getFilterMask(colors, categoryValues);
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

        const mask = this._getSelectionMask();
        const values = this._options.calculateColors(colors);

        let finalValues = values;
        if (mask.some((v) => !v)) {
            // We have some filtered out points.
            // If values are numbers, we need to convert to colors to mix with RGBA.
            if (values.length > 0 && typeof values[0] === 'number') {
                const numValues = values as number[];
                const range = this._getColorRange();
                const min = range.min;
                const max = range.max;

                finalValues = numValues.map((v, i) => {
                    if (!mask[i]) return '#d3d3d3';
                    return this._options.valueToColor(v, min, max);
                });
            } else {
                // Already strings (fixed color)
                finalValues = values.map((v, i) => {
                    if (!mask[i]) return '#d3d3d3';
                    return v;
                });
            }
        }

        // LOD: Apply filter to main trace values
        const mainValues =
            trace === 0 || trace === undefined ? this._applyLOD(finalValues) : finalValues;

        const selected = [];
        for (const data of this._selected.values()) {
            selected.push(data.color);
        }

        let dummy = [0.5, 0.5];
        if (this._options.hasColors()) {
            const range = this._getColorRange();
            dummy = [range.min, range.max];
        }

        return this._selectTrace<Array<string | number>>(mainValues, selected, dummy, trace);
    }

    /**
     * Get the numeric color values to use with the given plotly `trace`, or all of
     * them if `trace === undefined`. Used for customdata.
     */
    private _colorValues(trace?: number): Array<Array<number | string>> {
        if (this._options.hasColors()) {
            const prop = this._property(this._options.color.property.value);
            if (prop.string) {
                const interner = prop.string;
                const values = prop.values.map((v) => interner.string(v));

                const mainValues =
                    trace === 0 || trace === undefined ? this._applyLOD(values) : values;
                const selected = [];
                for (const data of this._selected.values()) {
                    const valIndex = prop.values[data.current];
                    selected.push(interner.string(valIndex));
                }
                const val0 = values[0] || '';
                const dummy = [val0, val0];
                return this._selectTrace(mainValues, selected, dummy, trace);
            } else {
                const rawValues = prop.values;
                const values = this._options.calculateColors(rawValues);

                const mainValues =
                    trace === 0 || trace === undefined ? this._applyLOD(values) : values;
                const selected = [];
                for (const data of this._selected.values()) {
                    selected.push(values[data.current]);
                }
                const range = this._getColorRange();
                const dummy = [range.min, range.max];
                return this._selectTrace(mainValues, selected, dummy, trace);
            }
        }

        const values = new Array<number>(
            this._property(this._options.x.property.value).values.length
        ).fill(NaN);
        const mainValues = trace === 0 || trace === undefined ? this._applyLOD(values) : values;
        return this._selectTrace(mainValues, [], [NaN, NaN], trace);
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
        let values = this._options.calculateSizes(sizes);
        const mask = this._getSelectionMask();
        if (mask.some((v) => !v)) {
            values = values.map((v, i) => (mask[i] ? v : v * 0.25));
        }

        // LOD: Apply filter to main trace values
        const mainValues = trace === 0 || trace === undefined ? this._applyLOD(values) : values;

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

        const dummy = [0, 0];

        return this._selectTrace<number | number[]>(
            mainValues as number | number[],
            selected,
            dummy,
            trace
        );
    }

    /**
     * Get the values to use as marker line color with the given plotly `trace`, or
     * all of them if `trace === undefined`.
     */
    private _lineColors(trace?: number): Array<string | string[]> {
        const n = this._property(this._options.x.property.value).values.length;
        // dafault to black outline
        let values = new Array(n).fill('black') as string[];

        const mask = this._getSelectionMask();
        // unselected points have transparent outline
        if (mask.some((v) => !v)) {
            values = values.map((v, i) => (mask[i] ? v : 'rgba(0,0,0,0)'));
        }

        // LOD: Apply filter to main trace values
        const mainValues = trace === 0 || trace === undefined ? this._applyLOD(values) : values;

        // Assume trace 1 (selection) line color is black.
        const selected = new Array<string>(this._selected.size).fill('black');

        // Dummy trace for colorbar is invisible
        const dummy = ['rgba(0,0,0,0)', 'rgba(0,0,0,0)'];

        return this._selectTrace<string | string[]>(
            mainValues as string | string[],
            selected,
            dummy,
            trace
        );
    }

    /**
     * Get the values to use as marker symbol with the given plotly `trace`, or
     * all of them if `trace === undefined`.
     */
    private _symbols(trace?: number): Array<string | string[] | number[]> {
        if (this._options.symbol.value === '') {
            // default to 0 (i.e. circles)
            return this._selectTrace<string | string[]>('circle', 'circle', 'circle', trace);
        }

        const property = this._property(this._options.symbol.value);
        const symbols = this._options.getSymbols(property);
        // LOD: Apply filter to main trace values
        const mainValues =
            trace === 0 || trace === undefined
                ? this._applyLOD(symbols as (string | number)[])
                : symbols;

        const selected = [];
        for (const data of this._selected.values()) {
            selected.push(symbols[data.current]);
        }

        const sym0 = symbols[0] || 'circle';
        const dummy = [sym0, sym0];

        return this._selectTrace<typeof symbols>(
            mainValues as typeof symbols,
            selected as typeof symbols,
            dummy as typeof symbols,
            trace
        );
    }

    /** Should we show the legend for the various symbols used? */
    private _showlegend(): boolean[] {
        const result = [false, false, false];

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
        const result = ['', '', ''];

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
     * Select either main, selected, dummy or all of them depending on `trace`, and return
     * them in a mode usable with `Plotly.restyle`/{@link PropertiesMap._restyle}
     */
    private _selectTrace<T>(main: T, selected: T, dummy: T, trace?: number): T[] {
        if (trace === 0) {
            return [main];
        } else if (trace === 1) {
            return [selected];
        } else if (trace === 2) {
            return [dummy];
        } else if (trace === undefined) {
            return [main, selected, dummy];
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

        // Ensure 3D traces are prepared
        // Show selected environments markers is handled by _getTraces logic for 3D

        // _switch3D logic for markers:
        for (const data of this._selected.values()) {
            data.toggleVisible(false);
        }
        this._updateMarkers();

        // We do a full react to update config (modebar)
        // Settings (like markers) are applied in _getTraces/_getLayout
        // which rely on _options and internal state. _options.z is enabled above.

        void this._react(this._getTraces(), this._getLayout(), this._getConfig());
    }

    /** Switch current plot from 3D back to 2D */
    private _switch2D(): void {
        assert(!this._is3D());
        this._options.z.disable();

        // Show selected environments markers
        for (const data of this._selected.values()) {
            data.toggleVisible(true);
        }

        // We do a full react to update config (modebar)
        // _getTraces uses _sizes() and options to set marker line width and size
        // consistently in 2D and 3D

        void this._react(this._getTraces(), this._getLayout(), this._getConfig());
    }

    /**
     * Get the Plotly configuration based on current mode (2D/3D)
     */
    private _getConfig(): Config {
        // Clone default config
        const config = { ...DEFAULT_CONFIG };
        config.modeBarButtonsToRemove = [...DEFAULT_CONFIG.modeBarButtonsToRemove];
        config.modeBarButtonsToAdd = [...DEFAULT_CONFIG.modeBarButtonsToAdd];

        if (this._is3D()) {
            config.modeBarButtonsToRemove.push('resetCameraDefault3d');
            // eslint-disable-next-line @typescript-eslint/no-unsafe-assignment
            config.modeBarButtonsToAdd = [
                [
                    {
                        name: 'Reset View',
                        // eslint-disable-next-line @typescript-eslint/no-explicit-any, @typescript-eslint/no-unsafe-assignment, @typescript-eslint/no-unsafe-member-access
                        icon: (Plotly as any).Icons.home,
                        click: () => {
                            this._resetToGlobalView();
                        },
                    },
                    // eslint-disable-next-line @typescript-eslint/no-explicit-any
                ] as any,
                ...config.modeBarButtonsToAdd,
            ];
        }

        return config as unknown as Config;
    }

    /**
     * Helper to trigger a full update of the plot traces (main, selected, dummy)
     * when data, styling, or selection changes.
     */
    private _restyleFull(layoutUpdate?: Partial<Layout>): void {
        const fullUpdate: Record<string, unknown> = {
            x: this._coordinates(this._options.x),
            y: this._coordinates(this._options.y),
            z: this._coordinates(this._options.z),
            'marker.color': this._colors(),
            'marker.size': this._sizes(),
            'marker.symbol': this._symbols(),
            'marker.line.color': this._lineColors(),
            hovertemplate: this._options.hovertemplate(),
            customdata: this._colorValues(),
            'marker.opacity': this._options.color.opacity.value / 100,
            visible: this._selectTrace(true, true, this._options.hasColors()),
        };

        const layout = layoutUpdate || {};

        if (layoutUpdate || this._is3D()) {
            if (this._is3D() && !('scene.camera' in layout)) {
                // explicitely preserve the camera position when updating the plot
                // to prevent it from snapping back to the default position
                // see https://github.com/lab-cosmo/chemiscope/issues/310
                const currentLayout = this._plot._fullLayout;
                if (currentLayout.scene !== undefined) {
                    // @ts-expect-error scene is defined in the layout
                    layout['scene.camera'] = currentLayout.scene.camera;
                }
            }

            void Plotly.update(
                this._plot,
                fullUpdate as unknown as Data,
                layout as Layout,
                [0, 1, 2]
            );
        } else {
            // Update main (0), selected (1) and dummy (2) traces
            // Use Plotly.restyle directly to allow awaiting (fixing synchronization issues)
            // while keeping the _restyle wrapper synchronous for legacy calls.
            void Plotly.restyle(this._plot, fullUpdate as unknown as Data, [0, 1, 2]);
        }
    }

    /**
     * Computes the subset of points to display based on spatial grid binning (LOD).
     */
    private _computeLOD(bounds?: Bounds): void {
        // check if LOD is enabled
        if (!this._options.useLOD.value) {
            this._lodIndices = null;
            return;
        }

        const xProp = this._options.x.property.value;
        const yProp = this._options.y.property.value;
        const zProp = this._options.z.property.value;

        const xValues = this._property(xProp).values;

        // Check threshold
        if (xValues.length <= PropertiesMap.LOD_THRESHOLD) {
            this._lodIndices = null;
            return;
        }

        const yValues = this._property(yProp).values;
        const is3D = this._is3D() && zProp !== '';
        const zValues = is3D ? this._property(zProp).values : null;

        const lodSet = new Set<number>();

        // Coarse pass
        // Compute a sparser "global" grid of points for the full range of the dataset
        // to show "something" when we rotate, pan or zoom

        const lodIndices = computeLODIndices(
            xValues,
            yValues,
            zValues,
            undefined,
            PropertiesMap.LOD_THRESHOLD / 10
        );

        for (const id of lodIndices) {
            lodSet.add(id);
        }

        // Fine pass
        // Do a higher resolution subsampling for the points that are actually visible
        const fineIndices =
            is3D && zValues && this._options.camera.value && bounds
                ? computeScreenSpaceLOD(
                      xValues,
                      yValues,
                      zValues,
                      this._options.camera.value,
                      bounds,
                      PropertiesMap.LOD_THRESHOLD / 2
                  )
                : computeLODIndices(xValues, yValues, zValues, bounds, PropertiesMap.LOD_THRESHOLD);

        for (const id of fineIndices) {
            lodSet.add(id);
        }

        this._lodIndices = Array.from(lodSet).sort((a, b) => a - b);
    }

    /**
     * Helper to filter data arrays using the computed LOD indices.
     */
    private _applyLOD<T>(values: T[]): T[] {
        if (this._lodIndices === null) {
            return values;
        }

        const result = new Array<T>(this._lodIndices.length);
        for (let i = 0; i < this._lodIndices.length; i++) {
            result[i] = values[this._lodIndices[i]];
        }
        return result;
    }

    /**
     * Resets the view to global bounds and recomputes LOD on the full dataset.
     * Used for double-click and autoscale events.
     */
    private _resetToGlobalView() {
        // Skip if another LOD update is in progress
        if (this._lodBusy) {
            return;
        }

        this._lodBusy = true;

        // Reset settings to Auto (NaN)
        this._options.x.min.value = NaN;
        this._options.x.max.value = NaN;
        this._options.y.min.value = NaN;
        this._options.y.max.value = NaN;
        this._options.z.min.value = NaN;
        this._options.z.max.value = NaN;

        // 1. Force global LOD computation
        this._computeLOD();

        // 2. Update data traces first (re-render points with global LOD)
        // We do this BEFORE relayout so that 'autorange' calculates bounds
        // based on the full dataset, not the sliced one.
        //
        // 3. Prepare Layout Update
        const layoutUpdate: Record<string, unknown> = {};

        if (this._is3D()) {
            // In 3D, 'autorange: true' resets camera AND axes.
            layoutUpdate['scene.xaxis.autorange'] = true;
            layoutUpdate['scene.yaxis.autorange'] = true;
            layoutUpdate['scene.zaxis.autorange'] = true;
            layoutUpdate['scene.aspectratio'] = { x: 1, y: 1, z: 1 };
            layoutUpdate['scene.camera'] = {
                center: { x: 0, y: 0, z: 0 },
                eye: { x: 1.25, y: 1.25, z: 1.25 },
                projection: { type: 'orthographic' },
                up: { x: 0, y: 0, z: 1 },
            };
        } else {
            // In 2D, we trigger autorange on standard axes
            layoutUpdate['xaxis.autorange'] = true;
            layoutUpdate['yaxis.autorange'] = true;
        }

        // 4. Force the view reset
        this._restyleFull(layoutUpdate as unknown as Partial<Layout>);
        // Manually trigger marker update for 2D mode
        if (!this._is3D()) {
            this._updateMarkers();
        }

        // Release lock
        this._lodBusy = false;
        this._afterplot();
    }

    /**
     * Function used as callback to update the axis ranges in settings after
     * the user changes zoom or range on the plot
     */
    private _afterplot(): void {
        if (this._lodBusy) {
            return;
        }

        // Set camera
        if (this._is3D()) {
            const { camera, aspectratio } = this._plot._fullLayout.scene;
            this._options.camera.setValue(plotlyToCamera({ camera, aspectratio }), 'DOM');
        }

        const bounds = this._getBounds();

        this._updateAxisSettings(bounds);

        // Update markers for 2D mode
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
                    'marker.line.color': this._lineColors(1),
                    x: this._coordinates(this._options.x, 1),
                    y: this._coordinates(this._options.y, 1),
                    z: this._coordinates(this._options.z, 1),
                } as Data,
                1
            );
        } else {
            // LOD: Retrieve raw values for marker positioning, skipping _coordinates
            // because _coordinates applies LOD filtering which might hide the marker's index
            const allX = this._property(this._options.x.property.value).values;
            const allY = this._property(this._options.y.property.value).values;
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
    private _getBounds(): Bounds {
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

    /**
     * Recomputes LOD indices based on current bounds and updates the plot.
     * Called when the view changes (zoom/pan/rotate).
     */
    private _updateLOD(bounds: Bounds): void {
        // Early exit if another LOD update is already in progress
        if (this._lodBusy) {
            return;
        }

        this._lodBusy = true;

        this._computeLOD(bounds);

        this._restyleFull();
        this._lodBusy = false;
    }

    /** Sync axis min/max so the UI reflects the current visible ranges in the plot */
    private _updateAxisSettings(bounds: Bounds): void {
        this._options.x.min.value = bounds.x[0];
        this._options.x.max.value = bounds.x[1];
        this._options.y.min.value = bounds.y[0];
        this._options.y.max.value = bounds.y[1];

        if (bounds.z !== undefined) {
            this._options.z.min.value = bounds.z[0];
            this._options.z.max.value = bounds.z[1];
        }
    }
}

/** Extract the data associated with the first `path` element in an SVG string */
function extractSvgPath(svg: string) {
    const doc = document.createElement('div');
    doc.innerHTML = svg;
    return doc.getElementsByTagName('path')[0].getAttribute('d');
}
