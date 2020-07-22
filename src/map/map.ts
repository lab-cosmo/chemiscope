/**
 * @packageDocumentation
 * @module map
 */

import assert from 'assert';

import Plotly from './plotly/plotly-scatter';
import {Config, Data, Layout, PlotlyScatterElement} from './plotly/plotly-scatter';

import {Property} from '../dataset';

import {EnvironmentIndexer, Indexes} from '../indexer';
import {OptionModificationOrigin, SavedSettings} from '../options';
import {GUID, PositioningCallback} from '../utils';
import {enumerate, getByID, getFirstKey, sendWarning} from '../utils';

import {MapData, NumericProperty} from './data';
import {MarkerData} from './markers';
import {AxisOptions, MapOptions} from './options';

import {COLOR_MAPS} from './colorscales';

const DEFAULT_LAYOUT = {
    // coloraxis is used for the markers
    coloraxis: {
        cmax: 0,
        cmin: 0,
        colorbar: {
            len: 1,
            thickness: 20,
            title: {
                text: '' as (undefined | string),
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
            title: '' as (undefined | string),
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
    ],
};

// in 3D mode, only strings are supported for 'marker.symbol', and only very few
// of them. See https://github.com/plotly/plotly.js/issues/4205 as the plotly
// issue tracking more symbols in 3D mode.
const POSSIBLE_SYMBOLS_IN_3D = ['circle', 'square', 'diamond', 'cross', 'x'];

function get3DSymbol(i: number): string {
    return POSSIBLE_SYMBOLS_IN_3D[i % POSSIBLE_SYMBOLS_IN_3D.length];
}

// get the max/min of an array. Math.min(...array) fails with very large arrays
function arrayMaxMin(values: number[]): {max: number, min: number} {
    let max = Number.NEGATIVE_INFINITY;
    let min = Number.POSITIVE_INFINITY;
    for (const value of values) {
        if (value > max) {
            max = value;
        }
        if (value < min) {
            min = value;
        }
    }
    assert(isFinite(min) && isFinite(max));
    return {max, min};
}

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
     * @param id         HTML id of the DOM element where the map should live
     * @param indexer    [[EnvironmentIndexer]] used to translate indexes from
     *                   environments index to structure/atom indexes
     * @param properties properties to be displayed
     */
    constructor(
        config: { id: string, settings: SavedSettings },
        indexer: EnvironmentIndexer,
        properties: { [name: string]: Property },
    ) {
        this._indexer = indexer;
        this.onselect = () => {};
        this.activeChanged = () => {};
        this._selected = new Map();

        this._root = getByID(config.id);
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
            config.settings,
        );
        this._colorReset = getByID<HTMLButtonElement>('chsp-color-reset');

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
    public select(indexes: Indexes) {
        if (this._active === undefined) {
            throw Error('tries to update selected environment, but there is no active marker');
        }

        const data = this._selected.get(this._active);
        assert(data !== undefined);

        // Plotly.restyle fires the plotly_click event, so ensure we only run
        // the update once.
        // https://github.com/plotly/plotly.js/issues/1025
        if (data.current !== indexes.environment) {
            data.current = indexes.environment;
            // Sets the active marker on this map
            this._updateSelectedMarker(data);
        }
    }

    /**
     * Set the marker with given GUID as the active marker.
     *
     * @param guid the GUID of the new active viewer
     */
    public setActive(guid: GUID) {
        if (this._active !== undefined) {
            const oldData = this._selected.get(this._active);
            assert(oldData !== undefined);
            oldData.marker.classList.toggle('chsp-active-structure', false);
        }

        this._active = guid;
        const data = this._selected.get(this._active);
        assert(data !== undefined);
        data.marker.classList.toggle('chsp-active-structure', true);

        if (this._is3D()) {
            this._restyle({'marker.size': this._sizes(1)} as Data, 1);
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

        const marker = document.createElement('div');
        this._root.appendChild(marker);
        marker.classList.add('chsp-structure-marker');
        if (guid === this._active) {
            marker.classList.toggle('chsp-active-structure', true);
        }
        marker.id = `chsp-selected-${guid}`;
        marker.onclick = () => {
            this.setActive(guid);
            const activeData = this._selected.get(guid);
            assert(activeData !== undefined);
            const activeIndexes = this._indexer.from_environment(activeData.current);
            this.activeChanged(guid, activeIndexes);
        };
        marker.style.backgroundColor = color;

        const data = {
            color: color,
            current: indexes.environment,
            marker: marker,
        };
        this._selected.set(guid, data);

        if (this._is3D()) {
            data.marker.style.display = 'none';
        } else {
            this._updateSelectedMarker(data);
        }

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
        assert(data.marker.parentNode !== null);
        data.marker.parentNode.removeChild(data.marker);

        this._selected.delete(guid);

        if (this._is3D()) {
            // We have to update all markers in 3D mode since we can not
            // update just the one that we removed
            this._updateAll3DMarkers();
        }
    }

    /**
     * Apply saved settings to the map.
     */
    public applySettings(settings: SavedSettings) {
        this._options.applySettings(settings);
    }

    /**
     * Save the values of the current settings in a way that an be used with
     * [[applySettings]] or saved to JSON.
     */
    public saveSettings(): SavedSettings {
        return this._options.saveSettings();
    }

    /** Forward to Ploty.restyle */
    private _restyle(data: Partial<Data>, traces?: number | number[]) {
        Plotly.restyle(this._plot, data, traces).catch((e) => setTimeout(() => { throw e; }));
    }

    /** Forward to Ploty.relayout */
    private _relayout(layout: Partial<Layout>) {
        Plotly.relayout(this._plot, layout).catch((e) => setTimeout(() => { throw e; }));
    }

    /** Add all the required callback to the settings */
    private _connectSettings() {
        // ======= x axis settings
        this._options.x.property.onchange = () => {
            const values = this._xValues();
            this._restyle({ x: values }, [0, 1]);
            this._relayout({
                'scene.xaxis.title': this._options.x.property.value,
                'xaxis.title': this._options.x.property.value,
            } as unknown as Layout);
        };

        this._options.x.scale.onchange = () => {
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
            const values = this._yValues();
            this._restyle({ y:  values}, [0, 1]);
            this._relayout({
                'scene.yaxis.title': this._options.y.property.value,
                'yaxis.title': this._options.y.property.value,
            } as unknown as Layout);
        };

        this._options.y.scale.onchange = () => {
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
        this._options.z.property.onchange = () => {
            const was3D = (this._plot as any)._fullData[0].type === 'scatter3d';
            if (this._options.z.property.value === '') {
                if (was3D) {
                    this._switch2D();
                }
            } else {
                if (!was3D) {
                    this._switch3D();
                }
            }

            const values = this._zValues();
            this._restyle({ z: values } as Data, [0, 1]);
            this._relayout({
                'scene.zaxis.title': this._options.z.property.value,
            } as unknown as Layout);
        };

        this._options.z.scale.onchange = () => {
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
        if (this._options.color.property.value) {
            this._options.color.enable();
            this._colorReset.disabled = false;

            const values = this._colors(0)[0] as number[];
            const {min, max} = arrayMaxMin(values);

            this._options.color.min.value = min;
            this._options.color.max.value = max;
        } else {
            this._options.color.min.disable();
            this._colorReset.disabled = true;

            this._options.color.min.value = 0;
            this._options.color.max.value = 0;
        }

        this._options.color.property.onchange = () => {
            if (this._options.color.property.value !== '') {
                this._options.color.enable();
                this._colorReset.disabled = false;

                const values = this._colors(0)[0] as number[];
                const {min, max} = arrayMaxMin(values);

                this._options.color.min.value = min;
                this._options.color.max.value = max;

                this._relayout({
                    'coloraxis.colorbar.title.text': this._options.color.property.value,
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

            this._restyle({
                'hovertemplate': this._hovertemplate(),
                'marker.color': this._colors(0),
            } as Data, 0);
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
                // seems to do the trick. This is possiblely a Ploty bug, we
                // would need to investiguate a bit more.
                'coloraxis.colorscale': this._colorScale(),
            } as unknown as Layout);
        };
        this._options.color.min.onchange = colorRangeChange;
        this._options.color.max.onchange = colorRangeChange;

        this._colorReset.onclick = () => {
            const values = this._colors(0)[0] as number[];
            const {min, max} = arrayMaxMin(values);
            this._options.color.min.value = min;
            this._options.color.max.value = max;
            this._relayout({
                'coloraxis.cmax': max,
                'coloraxis.cmin': min,
                // same as above regarding update of the points color
                'coloraxis.colorscale': this._colorScale(),
            } as unknown as Layout);
        };

        // ======= color palette
        this._options.palette.onchange = () => {
            this._relayout({
                'coloraxis.colorscale': this._colorScale(),
            } as unknown as Layout);
        };

        // ======= markers symbols
        this._options.symbol.onchange = () => {
            this._restyle({ 'marker.symbol': this._symbols() }, [0, 1]);

            this._restyle({
                name: this._legendNames(),
                showlegend: this._showlegend(),
            } as unknown as Data);

            this._relayout({
                'coloraxis.colorbar.len': this._colorbarLen(),
            } as unknown as Layout);
        };

        // ======= markers size
        const sizeChange = () => {
            this._restyle({ 'marker.size': this._sizes(0) } as Data, 0);
        };

        this._options.size.property.onchange = sizeChange;
        this._options.size.factor.onchange = sizeChange;
    }

    /** Actually create the Plotly plot */
    private _createPlot() {
        this._plot.innerHTML = '';

        const colors = this._colors();
        const lineColors = this._lineColors();
        const sizes = this._sizes();
        const symbols = this._symbols();

        const x = this._xValues();
        const y = this._yValues();
        const z = this._zValues();

        const type = this._is3D() ? 'scatter3d' : 'scattergl';

        // The main trace, containing default data
        const main = {
            name: '',
            type: type,

            x: x[0],
            y: y[0],
            z: z[0],

            hovertemplate: this._hovertemplate(),
            marker: {
                color: colors[0],
                coloraxis: 'coloraxis',
                line: {
                    color: lineColors[0],
                    width: 1,
                },
                // prevent plolty from messing with opacity when doing bubble
                // style charts (different sizes for each point)
                opacity: 1,
                size: sizes[0],
                sizemode: 'area',
                symbol: symbols[0],
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

            x: [NaN],
            y: [NaN],
            z: [NaN],

            hoverinfo: 'none',
            marker: {
                color: colors[1],
                line: {
                    color: lineColors[1],
                    width: 2,
                },
                opacity: 1,
                size: sizes[1],
                sizemode: 'area',
            },
            mode: 'markers',
            showlegend: false,
        };

        const traces = [main as Data, selected as Data];

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
        const layout = JSON.parse(JSON.stringify(DEFAULT_LAYOUT));
        // and set values speific to the displayed dataset
        layout.xaxis.title = this._options.x.property.value;
        layout.yaxis.title = this._options.y.property.value;
        layout.xaxis.type = this._options.x.scale.value;
        layout.yaxis.type = this._options.y.scale.value;
        layout.scene.xaxis.title = this._options.x.property.value;
        layout.scene.yaxis.title = this._options.y.property.value;
        layout.scene.zaxis.title = this._options.z.property.value;
        layout.coloraxis.colorscale = this._colorScale();
        layout.coloraxis.cmin = this._options.color.min.value;
        layout.coloraxis.cmax = this._options.color.max.value;
        layout.coloraxis.colorbar.title.text = this._options.color.property.value;
        layout.coloraxis.colorbar.len = this._colorbarLen();

        // Create an empty plot and fill it below
        Plotly.newPlot(this._plot, traces, layout as Partial<Layout>, DEFAULT_CONFIG as Config)
            .catch((e) => setTimeout(() => { throw e; }));

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
        this._updateAllMarkers();
    }

    /** Get the property with the given name */
    private _property(name: string): NumericProperty {
        const result = this._data[this._indexer.mode][name];
        if (result === undefined) {
            throw Error(`unknown property '${name}' requested in map`);
        }
        return result;
    }

    /** Get the plotly hovertemplate depending on `this._current.color` */
    private _hovertemplate(): string {
        if (this._hasColors()) {
            return this._options.color.property.value + ': %{marker.color:.2f}<extra></extra>';
        } else {
            return '%{x:.2f}, %{y:.2f}<extra></extra>';
        }
    }

    /**
     * Get the values to use for the x axis with the given plotly `trace`,
     * or all of them if `trace === undefined`
     */
    private _xValues(trace?: number): number[][] {
        const values = this._property(this._options.x.property.value).values;
        const selected = [];
        for (const marker of this._selected.values()) {
            if (this._is3D()) {
                selected.push(values[marker.current]);
            } else {
                selected.push(NaN);
            }
        }
        return this._selectTrace(values, selected, trace);
    }

    /**
     * Get the values to use for the y axis with the given plotly `trace`,
     * or all of them if `trace === undefined`
     */
    private _yValues(trace?: number): number[][] {
        const values = this._property(this._options.y.property.value).values;
        const selected = [];
        for (const marker of this._selected.values()) {
            if (this._is3D()) {
                selected.push(values[marker.current]);
            } else {
                selected.push(NaN);
            }
        }

        return this._selectTrace(values, selected, trace);
    }

    /**
     * Get the values to use for the z axis with the given plotly `trace`,
     * or all of them if `trace === undefined`
     */
    private _zValues(trace?: number): Array<undefined | number[]> {
        if (!this._is3D()) {
            return this._selectTrace(undefined, undefined, trace);
        }

        const values = this._property(this._options.z.property.value).values;
        const selected = [];
        for (const marker of this._selected.values()) {
            selected.push(values[marker.current]);
        }
        return this._selectTrace(values, selected, trace);
    }

    /**
     * Get the color values to use with the given plotly `trace`, or all of
     * them if `trace === undefined`
     */
    private _colors(trace?: number): Array<string | string[] | number | number[]> {
        let values;
        if (this._hasColors()) {
            values = this._property(this._options.color.property.value).values;
        } else {
            values = 0.5;
        }

        const selected = [];
        for (const data of this._selected.values()) {
            selected.push(data.color);
        }

        return this._selectTrace<string | string[] | number | number[]>(values, selected, trace);
    }

    /**
     * Get the **line** color values to use with the given plotly `trace`, or
     * all of them if `trace === undefined`
     */
    private _lineColors(trace?: number): string[] {
        if (this._is3D()) {
            return this._selectTrace<string>('black', 'black', trace);
        } else {
            return this._selectTrace<string>('rgba(1, 1, 1, 0.3)', 'black', trace);
        }
    }

    /** Get the colorscale to use for markers in the main plotly trace */
    private _colorScale(): Plotly.ColorScale {
        return COLOR_MAPS[this._options.palette.value];
    }

    /**
     * Get the values to use as marker size with the given plotly `trace`, or
     * all of them if `trace === undefined`.
     */
    private _sizes(trace?: number): Array<number | number[]> {
        // Transform the linear value from the slider into a logarithmic scale
        const logSlider = (value: number) => {
            const min_slider = 1;
            const max_slider = 100;

            // go from 1/6th of the size to 6 time the size
            const min_value = Math.log(1.0 / 6.0);
            const max_value = Math.log(6.0);

            const tmp = (max_value - min_value) / (max_slider - min_slider);
            return Math.exp(min_value + tmp * (value - min_slider));
        };

        const userFactor = logSlider(this._options.size.factor.value);

        let values;
        if (this._options.size.property.value !== '') {
            const sizes = this._property(this._options.size.property.value).values;
            const {min, max} = arrayMaxMin(sizes);
            const defaultSize = this._is3D() ? 2000 : 150;
            values = sizes.map((v: number) => {
                // normalize between 0 and 1, then scale by the user provided value
                const scaled = userFactor * (v + 0.05 - min) / (max - min);
                // since we are using scalemode: 'area', square the scaled value
                return defaultSize * scaled * scaled;
            });
        } else {
            // we need to use an array instead of a single value because of
            // https://github.com/plotly/plotly.js/issues/2735
            values = new Array(this._indexer.environmentsCount());
            if (this._is3D()) {
                values.fill(500 * userFactor);
            } else {
                values.fill(50 * userFactor);
            }
        }

        const selected = [];
        if (this._is3D()) {
            for (const guid of this._selected.keys()) {
                if (guid === this._active) {
                    selected.push(4000);
                } else {
                    selected.push(2000);
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

        const values = this._property(this._options.symbol.value).values;
        if (this._is3D()) {
            // If we need more symbols than available, we'll send a warning
            // and repeat existing ones
            if (this._symbolsCount() > POSSIBLE_SYMBOLS_IN_3D.length) {
              sendWarning(`${this._symbolsCount()} symbols are required, but we only have ${POSSIBLE_SYMBOLS_IN_3D.length}. Some symbols will be repeated`);
            }
            const symbols = values.map(get3DSymbol);
            const selected = [];
            for (const data of this._selected.values()) {
                selected.push(symbols[data.current]);
            }
            return this._selectTrace<string[]>(symbols, selected, trace);
        } else {
            // in 2D mode, use automatic assignment of symbols from numeric
            // values
            const selected = [];
            for (const data of this._selected.values()) {
                selected.push(values[data.current]);
            }
            return this._selectTrace<number[]>(values, selected, trace);
        }
    }

    /** Should we show the legend for the various symbols used? */
    private _showlegend(): boolean[] {
        const result = [false, false];

        if (this._options.symbol.value !== '') {
            for (let i = 0; i < this._symbolsCount()!; i++) {
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
            const names = this._property(this._options.symbol.value).string!.strings();
            for (const name of names) {
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

    /** How many different symbols are being displayed */
    private _symbolsCount(): number {
        if (this._options.symbol.value !== '') {
            return this._property(this._options.symbol.value).string!.strings().length;
        } else {
            return 0;
        }
    }

    /** Get the length of the colorbar to accomodate for the legend */
    private _colorbarLen(): number {
        /// Heigh of a legend item in plot unit
        const LEGEND_ITEM_HEIGH = 0.045;
        return 1 - LEGEND_ITEM_HEIGH * this._symbolsCount();
    }

    /** Does the current plot use color values? */
    private _hasColors(): boolean {
        return this._options.color.property.value !== '';
    }

    /** Is the the current plot a 3D plot? */
    private _is3D(): boolean {
        return this._options !== undefined && this._options.z.property.value !== '';
    }

    /** Switch current plot from 2D to 3D */
    private _switch3D() {
        assert(this._is3D());
        this._options.z.enable();

        const symbols = this._symbols();
        for (let s = 0; s < this._data.maxSymbols; s++) {
            symbols.push([get3DSymbol(s)]);
        }

        // switch all traces to 3D mode
        this._restyle({
            'marker.symbol': symbols,
            'type': 'scatter3d',
        } as unknown as Data);

        for (const data of this._selected.values()) {
            data.marker.style.display = 'none';
            this._updateSelectedMarker(data);
        }

        // Change the data that vary between 2D and 3D mode
        this._restyle({
            // transparency messes with depth sorting in 3D mode, even with
            // line width set to 0 ¯\_(ツ)_/¯
            // https://github.com/plotly/plotly.js/issues/4111
            'marker.line.color': this._lineColors(),
            'marker.line.width': [1, 2],
            // size change from 2D to 3D
            'marker.size': this._sizes(),
        } as Data, [0, 1]);

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
    private _switch2D() {
        assert(!this._is3D());
        this._options.z.disable();

        const symbols = this._symbols();
        for (let sym = 0; sym < this._data.maxSymbols; sym++) {
            symbols.push([sym]);
        }

        // switch all traces to 2D mode
        this._restyle({
            'marker.symbol': symbols,
            'type': 'scattergl',
        } as unknown as Data);

        // show selected environments markers
        for (const data of this._selected.values()) {
            data.marker.style.display = 'block';
        }

        // Change the data that vary between 2D and 3D mode
        this._restyle({
            // transparency messes with depth sorting in 3D mode
            // https://github.com/plotly/plotly.js/issues/4111
            'marker.line.color': this._lineColors(),
            'marker.line.width': [1, 0],
            // size change from 2D to 3D
            'marker.size': this._sizes(),
        } as Data, [0, 1]);

        this._relayout({
            // change colorbar length to accomodate for symbols legend
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
    private _afterplot() {
        // HACK: this is not public, so it might break
        const layout = this._plot._fullLayout;
        if (this._is3D()) {
            this._options.x.min.value = layout.scene.xaxis.range[0];
            this._options.x.max.value = layout.scene.xaxis.range[1];

            this._options.y.min.value = layout.scene.yaxis.range[0];
            this._options.y.max.value = layout.scene.yaxis.range[1];

            this._options.z.min.value = layout.scene.zaxis.range[0];
            this._options.z.max.value = layout.scene.zaxis.range[1];
        } else {
            this._options.x.min.value = layout.xaxis.range[0];
            this._options.x.max.value = layout.xaxis.range[1];

            this._options.y.min.value = layout.yaxis.range[0];
            this._options.y.max.value = layout.yaxis.range[1];

            this._updateAllMarkers();
        }
    }

    /**
     * Update the position of the given marker.
     *
     * In 3D mode, the markers uses the second Plotly trace.
     * In 2D mode, these markers are HTML div styled as colored circles that
     * we manually move around, saving a call to `restyle`.
     *
     * @param data data of the marker to update
     */
    private _updateSelectedMarker(data: MarkerData): void {
          const selected = data.current;
          const marker = data.marker;

          if (this._is3D()) {
              // we have to update all symbols at the same time
              this._updateAll3DMarkers();
          } else {
              const xaxis = this._plot._fullLayout.xaxis;
              const yaxis = this._plot._fullLayout.yaxis;

              const computeX = (value: number) => xaxis.l2p(value) + xaxis._offset;
              const computeY = (value: number) => yaxis.l2p(value) + yaxis._offset;

              const x = computeX(this._xValues(0)[0][selected]);
              const y = computeY(this._yValues(0)[0][selected]);

              // hide the point if it is outside the plot, allow for up to 10px
              // overflow (for points just on the border)
              const minX = computeX(xaxis.range[0]) - 10;
              const maxX = computeX(xaxis.range[1]) + 10;
              const minY = computeY(yaxis.range[1]) - 10;
              const maxY = computeY(yaxis.range[0]) + 10;
              if (!isFinite(x) || !isFinite(y) || x < minX || x > maxX || y < minY || y > maxY) {
                  marker.style.display = 'none';
              } else {
                  marker.style.display = 'block';
              }

              marker.style.top = `${y}px`;
              const plotWidth = this._plot.getBoundingClientRect().width;
              marker.style.right = `${plotWidth - x}px`;
          }
    }

    /**
     * Update the position, color & size of all markers in the map
     */
    private _updateAllMarkers(): void {
        if (this._is3D()) {
            this._updateAll3DMarkers();
        } else {
            for (const data of this._selected.values()) {
                this._updateSelectedMarker(data);
            }
        }
    }

    /**
     * Update the position, size & color for all markers in 3D mode
     */
    private _updateAll3DMarkers(): void {
        this._restyle({
            'marker.color': this._colors(1),
            'marker.size': this._sizes(1),
            'marker.symbol': this._symbols(1),
            'x': this._xValues(1),
            'y': this._yValues(1),
            'z': this._zValues(1),
        } as Data, 1);
    }
}
