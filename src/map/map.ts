/**
 * @packageDocumentation
 * @module map
 */

import assert from 'assert';

import Plotly from './plotly/plotly-scatter';
import { Config, Data, Layout, PlotlyScatterElement } from './plotly/plotly-scatter';

import { Property } from '../dataset';

import { EnvironmentIndexer, Indexes } from '../indexer';
import { OptionModificationOrigin, SavedSettings } from '../options';
import { GUID, PositioningCallback, arrayMaxMin, sendWarning } from '../utils';
import { enumerate, getByID, getElement, getFirstKey } from '../utils';

import { MapData, NumericProperty } from './data';
import { MarkerData } from './marker';
import { AxisOptions, MapOptions, get3DSymbol } from './options';

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
                text: '' as undefined | string,
            },
            y: 0,
            yanchor: 'bottom',
        },
        colorscale: [] as Plotly.ColorScale,
        showscale: true,
    },
    coloraxis2: {
        cmax: 0,
        cmin: 0,
        colorbar: {
            len: 1,
            thickness: 20,
            title: {
                text: '' as undefined | string,
            },
            y: 0,
            yanchor: 'bottom',
        },
        colorscale: [] as Plotly.ColorScale,
        showscale: false,
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
        range: undefined,
        title: '',
        type: 'linear',
        zeroline: false,
    },
    yaxis: {
        range: undefined,
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
                    Plotly.downloadImage(gd, {
                        filename: 'chemiscope-map',
                        format: 'png',
                        width: Math.max(gd._fullLayout.width, 600),
                        height: Math.max(gd._fullLayout.width, 600),
                    }).catch((e) =>
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
                    }).catch((e) =>
                        setTimeout(() => {
                            throw e;
                        })
                    );
                },
            },
        ],
    ],
};

const MAIN_TRACE: number = 0;
const SELECTED_TRACE: number = 1;
const BACKGROUND_TRACE: number = 2;
const ALL_TRACES: number[] = [MAIN_TRACE, SELECTED_TRACE, BACKGROUND_TRACE];
/**
 * The [[PropertiesMap]] class displays a 2D or 3D map (scatter plot) of
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

    /// HTML root holding the full plot
    private _root: HTMLElement;
    /// Plotly plot
    private _plot!: PlotlyScatterElement;
    /// All known properties
    private _data: MapData;

    /// Separation of Data into Foreground / Background
    private _is_foreground: Array<boolean>;

    /// GUID of the currently selected point
    private _active?: GUID;
    /// Map associating currently selected markers GUID to additional data
    private _selected: Map<GUID, MarkerData>;

    /// environment indexer
    private _indexer: EnvironmentIndexer;
    /// Settings of the map
    private _options: MapOptions;
    /// Button used to reset the range of color axis
    private _colorReset: HTMLButtonElement;

    /**
     * Create a new [[PropertiesMap]] inside the DOM element with the given HTML
     * `id`
     *
     * @param element   HTML element or string 'id' of the element where
     *                   the map should live
     * @param indexer    [[EnvironmentIndexer]] used to translate indexes from
     *                   environments index to structure/atom indexes
     * @param properties properties to be displayed
     */
    constructor(
        config: { element: string | HTMLElement; settings: SavedSettings },
        indexer: EnvironmentIndexer,
        properties: { [name: string]: Property }
    ) {
        this._indexer = indexer;
        this.onselect = () => {};
        this.activeChanged = () => {};
        this._selected = new Map<GUID, MarkerData>();

        this._root = getElement(config.element);

        if (this._root.style.position === '') {
            this._root.style.position = 'relative';
        }

        this._plot = document.createElement('div') as unknown as PlotlyScatterElement;
        this._plot.style.width = '100%';
        this._plot.style.height = '100%';
        this._root.appendChild(this._plot);

        this._data = new MapData(properties);

        this._options = new MapOptions(
            this._root,
            this._data[this._indexer.mode],
            (rect) => this.positionSettingsModal(rect),
            config.settings
        );
        this._colorReset = getByID<HTMLButtonElement>('chsp-color-reset');

        this._is_foreground = Array<boolean>(this._data.length(this._indexer.mode)).fill(true);
        this._connectSettings();

        // By default, position the modal for settings on top of the plot,
        // centered horizontally
        this.positionSettingsModal = (rect: DOMRect) => {
            const rootRect = this._root.getBoundingClientRect();
            return {
                left: rootRect.left + rootRect.width / 2 - rect.width / 2,
                top: rootRect.top + 20,
            };
        };

        this._createPlot();
    }

    /**
     * Remove all HTML added by this [[PropertiesMap]] in the current document
     */
    public remove(): void {
        this._root.innerHTML = '';
        this._options.remove();
    }

    /**
     * Change the environment indicated by the currently active marker to
     * the one with the given `indexes`
     */
    public select(indexes: Indexes): void {
        if (this._active === undefined) {
            throw Error('tries to update selected environment, but there is no active marker');
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
            this._restyle({ 'marker.size': this._sizes(SELECTED_TRACE) } as Data, SELECTED_TRACE);
        }
    }

    /**
     * Add a new marker to the map. The new marker is set as the active one
     *
     * @param guid GUID of the new marker
     * @param indexes indexes of the environment that the new marker should show
     */
    public addMarker(guid: GUID, color: string, indexes: Indexes): void {
        assert(!this._selected.has(guid));

        const data = new MarkerData(guid, color, indexes.environment);
        this._root.appendChild(data.marker);
        data.marker.onclick = () => {
            this.setActive(guid);
            this.activeChanged(guid, this._indexer.from_environment(data.current));
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
    public applySettings(settings: SavedSettings): void {
        this._options.applySettings(settings);
    }

    /**
     * Save the values of the current settings in a way that an be used with
     * [[applySettings]] or saved to JSON.
     */
    public saveSettings(): SavedSettings {
        return this._options.saveSettings();
    }

    /** Forward to Plotly.restyle */
    private _restyle(data: Partial<Data>, traces?: number | number[]) {
        Plotly.restyle(this._plot, data, traces).catch((e) =>
            setTimeout(() => {
                throw e;
            })
        );
    }

    /** Forward to Plotly.relayout */
    private _relayout(layout: Partial<Layout>) {
        Plotly.relayout(this._plot, layout).catch((e) =>
            setTimeout(() => {
                throw e;
            })
        );
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
                sendWarning(
                    'This property contains negative values. Note that taking the log will discard them.'
                );
            }
        };

        // ======= x axis settings
        this._options.x.property.onchange = () => {
            negativeLogWarning(this._options.x);
            const values = this._coordinates(this._options.x) as number[][];
            this._restyle({ x: values }, ALL_TRACES);
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
            this.setScaleStep(this._getBounds().x, 'x');
        };

        this._options.x.scale.onchange = () => {
            negativeLogWarning(this._options.x);
            this._options.setLogLabel(this._options.x, 'x');
            if (this._is3D()) {
                this._relayout({
                    'scene.xaxis.type': this._options.x.scale.value,
                } as unknown as Layout);
            } else {
                this._relayout({ 'xaxis.type': this._options.x.scale.value as Plotly.AxisType });
            }
        };

        // function creating a function to be used as onchange callback
        // for <axis>.min and <axis>.max
        const rangeChange = (name: string, axis: AxisOptions) => {
            return (_: number, origin: OptionModificationOrigin) => {
                if (origin === 'JS') {
                    // prevent recursion: this function calls relayout, which then
                    // calls _afterplot, which reset the min/max values.
                    return;
                }
                const min = axis.min.value;
                const max = axis.max.value;

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

        this._options.x.min.onchange = rangeChange('xaxis', this._options.x);
        this._options.x.max.onchange = rangeChange('xaxis', this._options.x);

        // ======= y axis settings
        this._options.y.property.onchange = () => {
            negativeLogWarning(this._options.y);
            const values = this._coordinates(this._options.y) as number[][];
            this._restyle({ y: values }, ALL_TRACES);
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
            this.setScaleStep(this._getBounds().y, 'y');
        };

        this._options.y.scale.onchange = () => {
            negativeLogWarning(this._options.y);
            this._options.setLogLabel(this._options.y, 'y');
            if (this._is3D()) {
                this._relayout({
                    'scene.yaxis.type': this._options.y.scale.value,
                } as unknown as Layout);
            } else {
                this._relayout({ 'yaxis.type': this._options.y.scale.value as Plotly.AxisType });
            }
        };

        this._options.y.min.onchange = rangeChange('yaxis', this._options.y);
        this._options.y.max.onchange = rangeChange('yaxis', this._options.y);

        // ======= z axis settings
        // setup initial state of the z axis settings
        if (this._options.z.property.value === '') {
            this._options.z.disable();
        } else {
            this._options.z.enable();
        }

        this._options.z.property.onchange = () => {
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
            this._restyle({ z: values } as Data, ALL_TRACES);
            this._relayout({
                'scene.zaxis.title': this._title(this._options.z.property.value),
                'scene.zaxis.autorange': true,
            } as unknown as Layout);
            if (this._is3D()) {
                this.setScaleStep(this._getBounds().z as number[], 'z');
            }
        };

        this._options.z.scale.onchange = () => {
            negativeLogWarning(this._options.z);
            this._options.setLogLabel(this._options.z, 'z');
            if (this._options.z.property.value !== '') {
                this._relayout({
                    'scene.zaxis.type': this._options.z.scale.value,
                } as unknown as Layout);
            }
        };

        this._options.z.min.onchange = rangeChange('zaxis', this._options.z);
        this._options.z.max.onchange = rangeChange('zaxis', this._options.z);

        // ======= color axis settings
        // setup initial state of the color settings
        if (this._options.color.property.value === '') {
            this._options.color.min.disable();
            this._colorReset.disabled = true;

            this._options.color.min.value = 0;
            this._options.color.max.value = 0;
        } else {
            this._options.color.enable();
            this._colorReset.disabled = false;

            const mainValues = this._colors(MAIN_TRACE)[0] as number[];
            const mainExt = arrayMaxMin(mainValues);

            const bgValues = this._colors(BACKGROUND_TRACE)[0] as number[];

            if (bgValues.length > 0) {
                const bgExt = arrayMaxMin(bgValues);
                this._options.color.min.value = Math.min(mainExt.min, bgExt.min);
                this._options.color.max.value = Math.max(mainExt.max, bgExt.max);
            } else {
                this._options.color.min.value = mainExt.min;
                this._options.color.max.value = mainExt.max;
            }
        }

        this._options.color.property.onchange = () => {
            if (this._options.color.property.value !== '') {
                this._options.color.enable();
                this._colorReset.disabled = false;

                const mainValues = this._colors(MAIN_TRACE)[0] as number[];
                const mainExt = arrayMaxMin(mainValues);

                const bgValues = this._colors(BACKGROUND_TRACE)[0] as number[];

                if (bgValues.length > 0) {
                    const bgExt = arrayMaxMin(bgValues);
                    this._options.color.min.value = Math.min(mainExt.min, bgExt.min);
                    this._options.color.max.value = Math.max(mainExt.max, bgExt.max);
                } else {
                    this._options.color.min.value = mainExt.min;
                    this._options.color.max.value = mainExt.max;
                }
                this._relayout({
                    'coloraxis.colorbar.title.text': this._title(
                        this._options.color.property.value
                    ),
                    'coloraxis.showscale': true,
                } as unknown as Layout);
            } else {
                this._options.color.disable();
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
                    'marker.color': this._colors(MAIN_TRACE),
                } as Data,
                MAIN_TRACE
            );
            this._restyle(
                {
                    hovertemplate: this._options.hovertemplate(),
                    'marker.color': this._colors(BACKGROUND_TRACE),
                } as Data,
                BACKGROUND_TRACE
            );
        };

        const colorRangeChange = () => {
            const min = this._options.color.min.value;
            const max = this._options.color.max.value;
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
            this._relayout({
                'coloraxis2.cmax': max,
                'coloraxis2.cmin': min,
                // looks like changing only 'coloraxis.cmax'/'coloraxis.cmin' do
                // not update the color of the points (although it does change
                // the colorbar). Asking for an update of 'coloraxis.colorscale'
                // seems to do the trick. This is possiblely a Ploty bug, we
                // would need to investiguate a bit more.
                'coloraxis2.colorscale': this._options.colorScale(),
            } as unknown as Layout);
        };
        this._options.color.min.onchange = colorRangeChange;
        this._options.color.max.onchange = colorRangeChange;

        this._colorReset.onclick = () => {
            const mainValues = this._colors(MAIN_TRACE)[0] as number[];
            const mainExt = arrayMaxMin(mainValues);

            const bgValues = this._colors(BACKGROUND_TRACE)[0] as number[];
            const bgExt = arrayMaxMin(bgValues);

            const min = Math.min(mainExt.min, bgExt.min);
            const max = Math.max(mainExt.max, bgExt.max);

            this._options.color.min.value = min;
            this._options.color.max.value = max;
            this._relayout({
                'coloraxis.cmax': max,
                'coloraxis.cmin': min,
                // same as above regarding update of the points color
                'coloraxis.colorscale': this._options.colorScale(),
            } as unknown as Layout);
            this._relayout({
                'coloraxis2.cmax': max,
                'coloraxis2.cmin': min,
                // same as above regarding update of the points color
                'coloraxis2.colorscale': this._options.colorScale(),
            } as unknown as Layout);
        };

        // ======= color palette
        this._options.palette.onchange = () => {
            this._relayout({
                'coloraxis.colorscale': this._options.colorScale(),
            } as unknown as Layout);
            this._relayout({
                'coloraxis2.colorscale': this._options.colorScale(),
            } as unknown as Layout);
        };

        // ======= opacity
        this._options.opacity.mode.onchange = () => {
            if (this._options.opacity.mode.value === 'filter') {
                this._options.opacity.enableColors.enable();
                this._options.opacity.filter.property.enable();
                this._options.opacity.filter.cutoff.enable();
                this._options.opacity.filter.operator.enable();
            } else {
                this._options.opacity.enableColors.disable();
                this._options.opacity.filter.property.disable();
                this._options.opacity.filter.cutoff.disable();
                this._options.opacity.filter.operator.disable();
            }

            this._updateTrace();

            if (this._options.opacity.mode.value === 'constant') {
                this._options.opacity.minimum.disable();
            } else {
                this._options.opacity.minimum.enable();
                this._options.opacity.minimum.value = Math.max(
                    0.0,
                    Math.min(
                        this._options.opacity.minimum.value,
                        this._options.opacity.maximum.value - 0.5
                    )
                );
                this._restyle(
                    {
                        'marker.opacity': this._options.opacity.minimum.value,
                    } as Data,
                    BACKGROUND_TRACE
                );
            }
            this._restyle(
                {
                    'marker.color': this._colors(MAIN_TRACE),
                    'marker.size': this._sizes(MAIN_TRACE),
                    'marker.symbol': this._symbols(MAIN_TRACE),
                    'marker.line.color': this._lineColors(MAIN_TRACE),
                    x: this._coordinates(this._options.x, MAIN_TRACE),
                    y: this._coordinates(this._options.y, MAIN_TRACE),
                    z: this._coordinates(this._options.z, MAIN_TRACE),
                } as Data,
                0
            );
            this._restyle(
                {
                    'marker.color': this._colors(BACKGROUND_TRACE),
                    'marker.size': this._sizes(BACKGROUND_TRACE),
                    'marker.symbol': this._symbols(BACKGROUND_TRACE),
                    'marker.line.color': this._lineColors(BACKGROUND_TRACE),
                    x: this._coordinates(this._options.x, BACKGROUND_TRACE),
                    y: this._coordinates(this._options.y, BACKGROUND_TRACE),
                    z: this._coordinates(this._options.z, BACKGROUND_TRACE),
                } as Data,
                BACKGROUND_TRACE
            );
        };

        this._options.opacity.filter.property.onchange = () => {
            this._updateTrace();
            this._restyle(
                {
                    'marker.color': this._colors(MAIN_TRACE),
                    'marker.size': this._sizes(MAIN_TRACE),
                    'marker.symbol': this._symbols(MAIN_TRACE),
                    'marker.line.color': this._lineColors(MAIN_TRACE),
                    x: this._coordinates(this._options.x, MAIN_TRACE),
                    y: this._coordinates(this._options.y, MAIN_TRACE),
                    z: this._coordinates(this._options.z, MAIN_TRACE),
                } as Data,
                MAIN_TRACE
            );
            this._restyle(
                {
                    'marker.color': this._colors(BACKGROUND_TRACE),
                    'marker.size': this._sizes(BACKGROUND_TRACE),
                    'marker.symbol': this._symbols(BACKGROUND_TRACE),
                    'marker.line.color': this._lineColors(BACKGROUND_TRACE),
                    x: this._coordinates(this._options.x, BACKGROUND_TRACE),
                    y: this._coordinates(this._options.y, BACKGROUND_TRACE),
                    z: this._coordinates(this._options.z, BACKGROUND_TRACE),
                } as Data,
                BACKGROUND_TRACE
            );
        };

        this._options.opacity.filter.operator.onchange = () => {
            this._updateTrace();
            this._restyle(
                {
                    'marker.color': this._colors(MAIN_TRACE),
                    'marker.size': this._sizes(MAIN_TRACE),
                    'marker.symbol': this._symbols(MAIN_TRACE),
                    'marker.line.color': this._lineColors(MAIN_TRACE),
                    x: this._coordinates(this._options.x, MAIN_TRACE),
                    y: this._coordinates(this._options.y, MAIN_TRACE),
                    z: this._coordinates(this._options.z, MAIN_TRACE),
                } as Data,
                MAIN_TRACE
            );
            this._restyle(
                {
                    'marker.color': this._colors(BACKGROUND_TRACE),
                    'marker.size': this._sizes(BACKGROUND_TRACE),
                    'marker.symbol': this._symbols(BACKGROUND_TRACE),
                    'marker.line.color': this._lineColors(BACKGROUND_TRACE),
                    x: this._coordinates(this._options.x, BACKGROUND_TRACE),
                    y: this._coordinates(this._options.y, BACKGROUND_TRACE),
                    z: this._coordinates(this._options.z, BACKGROUND_TRACE),
                } as Data,
                BACKGROUND_TRACE
            );
        };

        this._options.opacity.filter.cutoff.onchange = () => {
            this._updateTrace();
            this._restyle(
                {
                    'marker.color': this._colors(MAIN_TRACE),
                    'marker.size': this._sizes(MAIN_TRACE),
                    'marker.symbol': this._symbols(MAIN_TRACE),
                    'marker.line.color': this._lineColors(MAIN_TRACE),
                    x: this._coordinates(this._options.x, MAIN_TRACE),
                    y: this._coordinates(this._options.y, MAIN_TRACE),
                    z: this._coordinates(this._options.z, MAIN_TRACE),
                } as Data,
                MAIN_TRACE
            );
            this._restyle(
                {
                    'marker.color': this._colors(BACKGROUND_TRACE),
                    'marker.size': this._sizes(BACKGROUND_TRACE),
                    'marker.symbol': this._symbols(BACKGROUND_TRACE),
                    'marker.line.color': this._lineColors(BACKGROUND_TRACE),
                    x: this._coordinates(this._options.x, BACKGROUND_TRACE),
                    y: this._coordinates(this._options.y, BACKGROUND_TRACE),
                    z: this._coordinates(this._options.z, BACKGROUND_TRACE),
                } as Data,
                BACKGROUND_TRACE
            );
        };

        this._options.opacity.minimum.onchange = () => {
            this._restyle(
                {
                    'marker.opacity': this._options.opacity.minimum.value,
                } as Data,
                BACKGROUND_TRACE
            );
        };
        this._options.opacity.maximum.onchange = () => {
            this._restyle(
                {
                    'marker.opacity': this._options.opacity.maximum.value,
                } as Data,
                MAIN_TRACE
            );
        };

        this._options.opacity.enableColors.onchange = () => {
            this._restyle(
                {
                    'marker.color': this._colors(BACKGROUND_TRACE),
                } as Data,
                BACKGROUND_TRACE
            );
        };

        // ======= markers symbols
        this._options.symbol.onchange = () => {
            this._restyle({ 'marker.symbol': this._symbols() }, ALL_TRACES);

            this._restyle({
                name: this._legendNames(),
                showlegend: this._showlegend(),
            } as unknown as Data);

            this._relayout({
                'coloraxis.colorbar.len': this._colorbarLen(),
            } as unknown as Layout);
        };

        // ======= markers size
        // setup initial state of the marker size settings
        if (this._options.size.mode.value === '') {
            this._options.size.property.disable();
            this._options.size.reverse.disable();
        } else {
            this._options.size.property.enable();
            this._options.size.reverse.enable();
        }

        this._options.size.property.onchange = () => {
            if (this._options.size.property.value !== '') {
                this._options.size.mode.enable();
                this._options.size.reverse.enable();
            } else {
                this._options.size.mode.disable();
                this._options.size.reverse.disable();
            }
            this._restyle({ 'marker.size': this._sizes(MAIN_TRACE) } as Data, MAIN_TRACE);
            this._restyle(
                { 'marker.size': this._sizes(BACKGROUND_TRACE) } as Data,
                BACKGROUND_TRACE
            );
        };

        this._options.size.factor.onchange = () => {
            this._restyle({ 'marker.size': this._sizes(MAIN_TRACE) } as Data, MAIN_TRACE);
            this._restyle(
                { 'marker.size': this._sizes(BACKGROUND_TRACE) } as Data,
                BACKGROUND_TRACE
            );
        };

        this._options.size.mode.onchange = () => {
            this._restyle({ 'marker.size': this._sizes(MAIN_TRACE) } as Data, MAIN_TRACE);
            this._restyle(
                { 'marker.size': this._sizes(BACKGROUND_TRACE) } as Data,
                BACKGROUND_TRACE
            );
        };

        this._options.size.reverse.onchange = () => {
            this._restyle({ 'marker.size': this._sizes(MAIN_TRACE) } as Data, MAIN_TRACE);
            this._restyle(
                { 'marker.size': this._sizes(BACKGROUND_TRACE) } as Data,
                BACKGROUND_TRACE
            );
        };
    }

    /** Actually create the Plotly plot */
    private _createPlot() {
        this._plot.innerHTML = '';

        const type = this._is3D() ? 'scatter3d' : 'scattergl';
        // The main trace, containing default data
        const main = {
            name: '',
            type: type,

            x: this._coordinates(this._options.x, MAIN_TRACE)[0],
            y: this._coordinates(this._options.y, MAIN_TRACE)[0],
            z: this._coordinates(this._options.z, MAIN_TRACE)[0],

            hovertemplate: this._options.hovertemplate(),
            marker: {
                color: this._colors(MAIN_TRACE)[0],
                coloraxis: 'coloraxis',
                line: {
                    color: this._lineColors(MAIN_TRACE)[0],
                    width: 1,
                },
                // prevent plolty from messing with opacity when doing bubble
                // style charts (different sizes for each point)
                opacity: 1,
                size: this._sizes(MAIN_TRACE)[0],
                sizemode: 'area',
                symbol: this._symbols(MAIN_TRACE)[0],
            },
            mode: 'markers',
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

        // The backgraound trace, containing transparent data
        const main2 = {
            name: '',
            type: type,

            x: [0],
            y: [0],
            z: [0],

            hovertemplate: this._options.hovertemplate(),
            marker: {
                color: [0],
                coloraxis: 'coloraxis2',
                line: {
                    color: [0],
                    width: 1,
                },
                // prevent plolty from messing with opacity when doing bubble
                // style charts (different sizes for each point)
                opacity: 0,
                size: [0],
                sizemode: 'area',
                symbol: [0],
            },
            mode: 'markers',
            showlegend: false,
        };
        const traces = [main as Data, selected as Data, main2 as Data];

        const legendNames = this._legendNames().slice(2);
        const showlegend = this._showlegend().slice(2);
        assert(legendNames.length === showlegend.length);
        const currentLength = legendNames.length;

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
        layout.coloraxis.colorbar.title.text = this._title(this._options.color.property.value);
        layout.coloraxis.colorbar.len = this._colorbarLen();
        layout.coloraxis2.colorscale = this._options.colorScale();
        layout.coloraxis2.cmin = this._options.color.min.value;
        layout.coloraxis2.cmax = this._options.color.max.value;

        // Create an empty plot and fill it below
        Plotly.newPlot(
            this._plot,
            traces,
            layout as Partial<Layout>,
            DEFAULT_CONFIG as unknown as Config
        ).catch((e) =>
            setTimeout(() => {
                throw e;
            })
        );
        this._plot.classList.add('chsp-map');

        this._plot.on('plotly_click', (event: Plotly.PlotMouseEvent) => {
            // don't update selected env on double click, since it is bound to
            // 'reset zoom level' in 2D mode.
            if (event.event && event.event.detail === 2) {
                return;
            }

            let environment = event.points[0].pointNumber;
            if (this._is3D() && event.points[0].data.name === 'selected') {
                // if someone has clicked on a selection marker, set to active
                // this is only used in 3D mode, since in 2D the HTML marker
                // directly deal with the click event
                for (const [i, [guid, data]] of enumerate(this._selected.entries())) {
                    if (event.points[0].pointNumber === i) {
                        environment = data.current;
                        if (this._active !== guid) {
                            this.setActive(guid);
                            this.activeChanged(guid, this._indexer.from_environment(data.current));
                        }
                        break;
                    }
                }
            }

            const indexes = this._indexer.from_environment(environment);

            this.select(indexes);
            this.onselect(indexes);
        });

        this._plot.on('plotly_afterplot', () => this._afterplot());
        this._updateMarkers();

        // set step of min/max select arrows based on the plot range
        const bounds = this._getBounds();
        this.setScaleStep(bounds.x, 'x');
        this.setScaleStep(bounds.y, 'y');
        if (bounds.z !== undefined) {
            this.setScaleStep(bounds.z, 'z');
        }
    }

    /** Get the property with the given name */
    private _property(name: string): NumericProperty {
        const result = this._data[this._indexer.mode][name];
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
        const [main, background] = this._filterForOpacity<number>(values);
        return this._selectTrace<number[]>(main, selected, background, trace);
    }

    private _title(name: string): string {
        if (name !== '') {
            const units = this._property(name).units;
            if (units !== undefined) {
                return name + `/${units}`;
            }
        }
        return name;
    }

    /**
     * Get the color values to use with the given plotly `trace`, or all of
     * them if `trace === undefined`
     */
    private _colors(trace?: number): Array<string | string[] | number | number[]> {
        let values;
        if (this._options.hasColors()) {
            values = this._property(this._options.color.property.value).values;
        } else {
            values = 0.5;
        }

        let main, background;
        if (!this._options.opacity.enableColors.value) {
            [main, background] = this._filterForOpacity<number>(values as number[]);
            const a = 0.3 * this._options.opacity.minimum.value;
            background = `rgba(${a}, ${a}, ${a}, ${a})`;
        } else {
            [main, background] = this._filterForOpacity<number>(values as number[]);
        }

        const selected = [];
        for (const data of this._selected.values()) {
            selected.push(data.color);
        }

        return this._selectTrace<string | string[] | number | number[]>(
            main,
            selected,
            background,
            trace
        );
    }

    private _filterForOpacity<T>(values: Array<T>): Array<T[]> {
        const foreground: Array<T> = [];
        const background: Array<T> = [];

        assert(values.length === this._is_foreground.length);
        for (let i = 0; i < this._is_foreground.length; i++) {
            if (this._is_foreground[i]) {
                foreground.push(values[i]);
            } else {
                background.push(values[i]);
            }
        }
        return [foreground, background];
    }

    private _updateTrace(): void {
        if (this._options.opacity.mode.value === 'filter') {
            const filter_property = this._property(
                this._options.opacity.filter.property.value
            ).values;
            const filter_cutoff = this._options.opacity.filter.cutoff.value;
            const op = this._options.opacity.filter.operator.value;
            for (let i = 0; i < filter_property.length; i++) {
                if (op === '=' && filter_property[i] !== filter_cutoff) {
                    this._is_foreground[i] = false;
                } else if (op === '<' && filter_property[i] >= filter_cutoff) {
                    this._is_foreground[i] = false;
                } else if (op === '>' && filter_property[i] <= filter_cutoff) {
                    this._is_foreground[i] = false;
                } else {
                    this._is_foreground[i] = true;
                }
            }
        } else {
            this._is_foreground = Array<boolean>(this._data.length(this._indexer.mode)).fill(true);
        }
    }

    /**
     * Get the **line** color values to use with the given plotly `trace`, or
     * all of them if `trace === undefined`
     */
    private _lineColors(trace?: number): Array<string | string[]> {
        if (this._is3D()) {
            return this._selectTrace<string>('black', 'black', 'rgba(0,0,0,0)', trace);
        } else {
            const defaultColor = `rgba(1, 1, 1, ${0.3 * this._options.opacity.maximum.value})`;
            const lighterColor = `rgba(1, 1, 1, ${0.3 * this._options.opacity.minimum.value})`;
            return this._selectTrace<string>(defaultColor, defaultColor, lighterColor, trace);
        }
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
            );
        }
        const values = this._options.calculateSizes(sizes);
        const [main, background] = this._filterForOpacity<number>(values);

        const selected = [];
        if (this._is3D()) {
            for (const guid of this._selected.keys()) {
                if (guid === this._active) {
                    selected.push(600);
                } else {
                    selected.push(300);
                }
            }
        }
        return this._selectTrace<number | number[]>(main, selected, background, trace);
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
        const [main, background] = this._filterForOpacity<string>(symbols as string[]);
        const selected = [];
        for (const data of this._selected.values()) {
            selected.push(symbols[data.current]);
        }
        return this._selectTrace<typeof symbols>(
            main,
            selected as typeof symbols,
            background,
            trace
        );
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
     * them in a mode usable with `Plotly.restyle`/[[PropertiesMap._restyle]]
     */
    private _selectTrace<T>(main: T, selected: T, background?: T, trace?: number): T[] {
        if (background === undefined) {
            if (Array.isArray(main)) {
                [main, background] = this._filterForOpacity<typeof main[0]>(main) as unknown as T[];
            } else {
                background = main;
            }
        }
        if (trace === MAIN_TRACE) {
            return [main];
        } else if (trace === SELECTED_TRACE) {
            return [selected];
        } else if (trace === BACKGROUND_TRACE) {
            return [background];
        } else if (trace === undefined) {
            return [main, selected, background];
        } else {
            throw Error('internal error: invalid trace number');
        }
    }

    /** Get the length of the colorbar to accommodate for the legend */
    private _colorbarLen(): number {
        /// Heigh of a legend item in plot unit
        const LEGEND_ITEM_HEIGH = 0.045;
        return 1 - LEGEND_ITEM_HEIGH * this._symbolsCount();
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
        this._restyle(
            {
                // transparency messes with depth sorting in 3D mode, even with
                // line width set to 0 ¯\_(ツ)_/¯
                // https://github.com/plotly/plotly.js/issues/4111
                'marker.line.color': this._lineColors(),
                'marker.line.width': [1, 2, 1],
                // size change from 2D to 3D
                'marker.size': this._sizes(),
                'marker.sizemode': 'area',
            } as Data,
            ALL_TRACES
        );

        if (
            this._options.opacity.minimum.value !== 0.0 ||
            this._options.opacity.maximum.value !== 1.0
        ) {
            sendWarning('In 3D all markers must be either transparent or opaque.');
        }
        this._options.opacity.minimum.disable();
        this._options.opacity.maximum.disable();
        this._options.opacity.enableColors.disable();
        this._restyle({ opacity: 0 } as Data, [BACKGROUND_TRACE]);

        this._relayout({
            // change colorbar length to accomodate for symbols legend
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

        if (
            this._options.opacity.minimum.value !== 0.0 ||
            this._options.opacity.maximum.value !== 1.0
        ) {
            sendWarning('Resetting original maximum and minimum opacities.');
        }
        this._options.opacity.minimum.enable();
        this._options.opacity.maximum.enable();
        this._options.opacity.enableColors.enable();
        this._restyle(
            { opacity: 1, 'marker.opacity': this._options.opacity.maximum.value } as Data,
            [MAIN_TRACE]
        );
        this._restyle(
            { opacity: 1, 'marker.opacity': this._options.opacity.minimum.value } as Data,
            [BACKGROUND_TRACE]
        );

        this._restyle(
            {
                // transparency messes with depth sorting in 3D mode
                // https://github.com/plotly/plotly.js/issues/4111
                'marker.line.color': this._lineColors(),
                'marker.line.width': [1, 0, 1],
                // size change from 2D to 3D
                'marker.size': this._sizes(),
            } as Data,
            ALL_TRACES
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

        this._options.x.min.value = bounds.x[0];
        this._options.x.max.value = bounds.x[1];
        this._options.y.min.value = bounds.y[0];
        this._options.y.max.value = bounds.y[1];
        if (bounds.z !== undefined) {
            this._options.z.min.value = bounds.z[0];
            this._options.z.max.value = bounds.z[1];
        }

        if (!this._is3D()) {
            this._updateMarkers();
        }
    }

    /**
     * Update the position, color & size of markers within the data array
     */
    private _updateMarkers(data: MarkerData[] = Array.from(this._selected.values())): void {
        if (this._is3D()) {
            data.forEach((d) => d.toggleVisible(false));
            this._restyle(
                {
                    'marker.color': this._colors(SELECTED_TRACE),
                    'marker.size': this._sizes(SELECTED_TRACE),
                    'marker.symbol': this._symbols(SELECTED_TRACE),
                    x: this._coordinates(this._options.x, SELECTED_TRACE),
                    y: this._coordinates(this._options.y, SELECTED_TRACE),
                    z: this._coordinates(this._options.z, SELECTED_TRACE),
                } as Data,
                SELECTED_TRACE
            );
        } else {
            const allX = this._coordinates(this._options.x, MAIN_TRACE) as number[][];
            const allY = this._coordinates(this._options.y, MAIN_TRACE) as number[][];
            const plotWidth = this._plot.getBoundingClientRect().width;

            for (const datum of data) {
                const rawX = allX[0][datum.current];
                const rawY = allY[0][datum.current];
                if (this._insidePlot(rawX, rawY)) {
                    const x = plotWidth - this._computeRSCoord(rawX, 'x');
                    const y = this._computeRSCoord(rawY, 'y');
                    datum.update(x, y);
                } else {
                    datum.toggleVisible();
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

    // Computes the real space coordinate of a value on the plot
    private _computeRSCoord(value: number, axisName: string): number {
        assert(axisName === 'x' || axisName === 'y');
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

    // Checks if a point is in the visible plot
    private _insidePlot(x: number, y: number, z?: number): boolean {
        const tolerance = 10;

        const bounds = this._getBounds();
        const isInsideRange = (value: number, range: [number, number], tolerance: number) => {
            return value > range[0] - tolerance && value < range[1] + tolerance;
        };

        let inside = isInsideRange(x, bounds.x, tolerance);
        inside = inside && isInsideRange(y, bounds.y, tolerance);

        if (z !== undefined) {
            assert(bounds.z !== undefined);
            inside = inside && isInsideRange(z, bounds.z, tolerance);
        }
        return inside;
    }

    /** Changes the step of the arrow buttons in min/max input based on dataset range*/
    private setScaleStep(axisBounds: number[], axisName: string): void {
        if (axisBounds !== undefined) {
            // round to 10 decimal places so it does not break in Firefox
            const step = Math.round(((axisBounds[1] - axisBounds[0]) / 20) * 10 ** 10) / 10 ** 10;
            const minElement = getByID<HTMLInputElement>(`chsp-${axisName}-min`);
            const maxElement = getByID<HTMLInputElement>(`chsp-${axisName}-max`);
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
