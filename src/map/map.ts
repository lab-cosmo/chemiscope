/**
 * @packageDocumentation
 * @module map
 */

import assert from "assert";

import {Property} from "../dataset";
import {makeDraggable, getByID, Indexes, EnvironmentIndexer} from "../utils";

import * as Plotly from "./plotly/plotly-scatter";
import {Config, Data, Layout, PlotlyScatterElement} from "./plotly/plotly-scatter";

import {COLOR_MAPS} from "./colorscales";
import {MapData, NumericProperty} from "./data";

import HTML_SETTINGS from "./settings.html";

const DEFAULT_LAYOUT = {
    hovermode: "closest",
    showlegend: true,
    legend: {
        y: 1,
        yanchor: "top",
        tracegroupgap: 5,
        itemclick: false,
        itemdoubleclick: false,
    },
    title: {
        text: "",
        y: 0.982,
        x: 0.08,
        font: {
            size: 25,
        },
    },
    margin: {
        l: 50,
        t: 50,
        r: 50,
        b: 50,
    },
    xaxis: {
        title: "",
        zeroline: false,
        type: "linear",
        range: undefined,
    },
    yaxis: {
        title: "",
        zeroline: false,
        type: "linear",
        range: undefined,
    },
    // coloraxis is used for the markers
    coloraxis: {
        cmin: 0,
        cmax: 0,
        colorscale: <Plotly.ColorScale> [],
        showscale: true,
        colorbar: {
            title: {
                text: <undefined | string> "",
            },
            thickness: 20,
            len: 1,
            y: 0,
            yanchor: "bottom",
        }
    },
    // coloraxis2 is used for the markers lines
    coloraxis2: {
        cmin: 0,
        cmax: 0,
        colorscale: <Plotly.ColorScale> [],
        showscale: false,
    },
    scene: {
        camera: {
            projection: {
                // broken, see https://github.com/plotly/plotly.js/issues/4514
                // type: "orthographic"
            }
        },
        xaxis: {
            title: "",
            showspikes: false,
        },
        yaxis: {
            title: "",
            showspikes: false,
        },
        zaxis: {
            title: <undefined | string> "",
            showspikes: false,
        },
    },
};

const DEFAULT_CONFIG = {
    displayModeBar: true,
    displaylogo: false,
    scrollZoom: true,
    modeBarButtonsToRemove: [
        "hoverClosestCartesian",
        "hoverCompareCartesian",
        "toggleSpikelines",
        "autoScale2d",
        "zoomIn2d",
        "zoomOut2d",
        "select2d",
        "lasso2d",
        "hoverClosest3d",
        "tableRotation",
        "resetCameraLastSave3d",
    ],
};

/** HTML element holding settings for a given axis (x, y, z, color) */
interface AxisSetting {
    /// Which values shoudl we use
    select: HTMLSelectElement;
    /// Which scale (linear/log) should we use
    scale: HTMLSelectElement;
    /// The minimal value for this axis
    min: HTMLInputElement;
    /// The maximal value for this axis
    max: HTMLInputElement;
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
    /// HTML root holding the full plot
    private _root: HTMLElement;
    /// Plotly plot
    private _plot!: PlotlyScatterElement;
    /// The dataset name
    private _name: string;
    /// All known properties
    private _data: MapData;
    /// Index of the currently selected point
    private _selected: number;
    /// environment indexer
    private _indexer: EnvironmentIndexer;
    /// Currently displayed data
    private _current!: {
        /// Name of the properties in `this._properties()` used for x values
        x: string,
        /// Name of the properties in `this._properties()` used for y values
        y: string,
        /// Name of the properties in `this._properties()` used for z values,
        /// `undefined` for 2D plots
        z?: string,
        /// Name of the colorscale to use
        colorscale: string,
        /// Name of the properties in `this._properties()` used for color
        /// values, `undefined` when using the default colors
        color?: string,
        /// Name of the properties in `this._properties()` used for size values,
        /// `undefined` when using the default sizes
        size?: string,
        /// Name of the properties in `this._properties()` used for symbols
        /// values, `undefined` when using the default symbols
        symbols?: string,
    }
    /// callback to get the initial positioning of the settings modal. The
    /// callback gets the current placement of the settings as a DOMRect, and
    /// should return top and left positions in pixels, used with
    /// `position: fixed`
    private _settingsPlacement!: (rect: DOMRect) => {top: number, left: number};
    /// Store the HTML elements used for settings
    private _settings: {
        x: AxisSetting;
        y: AxisSetting;
        z: AxisSetting;
        color: AxisSetting;
        colorReset: HTMLButtonElement;
        palette: HTMLSelectElement;
        symbol: HTMLSelectElement;
        size: {
            select: HTMLSelectElement;
            factor: HTMLInputElement;
        };
    };
    /// Marker indicating the position of the latest selected point in 2D mode
    /// Using such div is much faster than trying to restyle the full plot,
    /// especially with more than 100k points. In 3D mode, a separate trace is
    /// used instead.
    private _selectedMarker: HTMLElement;

    /** Callback fired when the plot is clicked and a new point is selected */
    public onselect: (indexes: Indexes) => void;

    /**
     * Create a new [[PropertiesMap]] inside the DOM element with the given HTML
     * `id`
     *
     * @param id         HTML id of the DOM element where the map should live
     * @param name       name of the dataset, displayed as the plot title
     * @param indexer    [[EnvironmentIndexer]] used to translate indexes from
     *                   environments index to structure/atom indexes
     * @param properties properties to be displayed
     */
    constructor(id: string, name: string, indexer: EnvironmentIndexer, properties: {[name: string]: Property}) {
        this._name = name;
        this._indexer = indexer;
        this._selected = 0;
        this.onselect = () => {};

        const root = document.getElementById(id);
        if (root === null) {
            throw Error(`could not find HTML element #${id}`)
        }
        this._root = root;
        this._root.style.position = 'relative';

        this._plot = document.createElement("div") as unknown as PlotlyScatterElement;
        this._plot.style.width = "100%";
        this._plot.style.height = "100%";
        this._plot.style.minHeight = "550px";
        this._root.appendChild(this._plot);

        this._selectedMarker = document.createElement("div");
        this._root.appendChild(this._selectedMarker);
        this._selectedMarker.classList.add("chsp-selected");

        this._data = new MapData(properties);
        this._setupDefaults();

        this._createSettings();
        this._settings = {
            x: {
                select: getByID<HTMLSelectElement>('chsp-x'),
                scale: getByID<HTMLSelectElement>('chsp-x-scale'),
                min: getByID<HTMLInputElement>('chsp-x-min'),
                max: getByID<HTMLInputElement>('chsp-x-max'),
            },
            y: {
                select: getByID<HTMLSelectElement>('chsp-y'),
                scale: getByID<HTMLSelectElement>('chsp-y-scale'),
                min: getByID<HTMLInputElement>('chsp-y-min'),
                max: getByID<HTMLInputElement>('chsp-y-max'),
            },
            z: {
                select: getByID<HTMLSelectElement>('chsp-z'),
                scale: getByID<HTMLSelectElement>('chsp-z-scale'),
                min: getByID<HTMLInputElement>('chsp-z-min'),
                max: getByID<HTMLInputElement>('chsp-z-max'),
            },
            color: {
                select: getByID<HTMLSelectElement>('chsp-color'),
                scale: document.createElement('select'),
                min: getByID<HTMLInputElement>('chsp-color-min'),
                max: getByID<HTMLInputElement>('chsp-color-max'),
            },
            colorReset: getByID<HTMLButtonElement>('chsp-color-reset'),
            palette: getByID<HTMLSelectElement>('chsp-palette'),
            symbol: getByID<HTMLSelectElement>('chsp-symbol'),
            size: {
                select: getByID<HTMLSelectElement>('chsp-size'),
                factor: getByID<HTMLInputElement>('chsp-size-factor'),
            }
        };

        this._connectSettings();
        this._setupSettings();

        this._createPlot();
    }

    /** Change the selected environment to the one with the given `indexes` */
    public select(indexes: Indexes) {
        // Plotly.restyle fires the plotly_click event, so ensure we only run
        // the update once.
        // https://github.com/plotly/plotly.js/issues/1025
        if (indexes.environment !== this._selected) {
            this._selected = indexes.environment;
            this._updateSelectedMarker();
        }
    }

    /**
     * Change the displayed dataset to a new one, without re-creating the
     * viewer itself.
     *
     * @param  name       name of the new dataset
     * @param  indexer    new indexer making the environment index to
     *                    structure/atom pair translation
     * @param  properties new properties to display
     */
    public changeDataset(name: string, indexer: EnvironmentIndexer, properties: {[name: string]: Property}) {
        if (this._is3D()) {
            this._current.z = undefined;
            this._switch2D();
        }

        this._name = name;
        this._indexer = indexer;
        this._selected = 0;
        this._data = new MapData(properties);
        this._setupDefaults();
        this._setupSettings();
        this._createPlot();

        this._relayout({ 'xaxis.autorange': true, 'yaxis.autorange': true });
    }

    /**
     * Use the given callback to compute the placement of the settings modal.
     * The callback gets the current placement of the settings as a DOMRect,
     * and should return top and left positions in pixels, used with `position:
     * fixed`. The callback is called once, the first time the settings are
     * opened.
     */
    public settingsPlacement(callback: (rect: DOMRect) => {top: number, left: number}) {
        this._settingsPlacement = callback;
    }

    /** Forward to Ploty.restyle */
    private _restyle(data: Partial<Data>, traces?: number | number[]) {
        Plotly.restyle(this._plot, data, traces).catch(e => setTimeout(() => {throw e;}));
    }

    /** Forward to Ploty.relayout */
    private _relayout(layout: Partial<Layout>) {
        Plotly.relayout(this._plot, layout).catch(e => setTimeout(() => {throw e;}));
    }

    /** setup the default values after loading a new dataset */
    private _setupDefaults() {
        const prop_names = Object.keys(this._properties());
        if (prop_names.length < 2) {
            throw Error("Sketchmap needs at least two properties to plot")
        }

        this._current = {
            x: prop_names[0],
            y: prop_names[1],
            colorscale: 'inferno',
        }

        if (prop_names.length > 2) {
            this._current.color = prop_names[2];
        }
    }

    /** Create the settings modal by adding HTML to the page */
    private _createSettings() {
        // use HTML5 template to generate a DOM object from an HTML string
        const template = document.createElement('template');
        template.innerHTML = `<button data-target="#chsp-settings"
                                      data-toggle="modal"
                                      class="btn btn-light btn-sm chsp-open-map-settings">
            <div class="chsp-hamburger"><div></div><div></div><div></div></div>
        </button>`;
        const openSettings = template.content.firstChild!;
        this._root.append(openSettings);

        // replace id to ensure they are unique even if we have mulitple viewers
        // on a single page
        template.innerHTML = HTML_SETTINGS;
        const modal = template.content.firstChild!;
        document.body.appendChild(modal);

        const modalDialog = modal.childNodes[1]! as HTMLElement
        if (!modalDialog.classList.contains('modal-dialog')) {
            throw Error("internal error: missing modal-dialog class")
        }
        // make the settings modal draggable
        makeDraggable(modalDialog, ".modal-header");

        // Position modal near the actual viewer
        openSettings.addEventListener('click', () => {
            // only set style once, on first open, and keep previous position
            // on next open to keep the 'draged-to' position
            if (modalDialog.getAttribute('data-initial-modal-positions-set') === null) {
                modalDialog.setAttribute('data-initial-modal-positions-set', "true");

                // display: block to ensure modalDialog.offsetWidth is non-zero
                (modalDialog.parentNode as HTMLElement).style.display = 'block';

                const {top, left} = this._settingsPlacement(modalDialog.getBoundingClientRect());

                // set width first, since setting position can influence it
                modalDialog.style.width = `${modalDialog.offsetWidth}px`;
                // unset margins when using position: fixed
                modalDialog.style.margin = '0';
                modalDialog.style.position = 'fixed';
                modalDialog.style.top = `${top}px`;
                modalDialog.style.left = `${left}px`;
            }
        })

        // By default, position the modal for settings on top of the plot,
        // centered horizontally
        this._settingsPlacement = (rect: DOMRect) => {
            const rootRect = this._root.getBoundingClientRect();
            return {
                top: rootRect.top + 20,
                left: rootRect.left + rootRect.width / 2 - rect.width / 2,
            }
        }
    }

    /** Add all the required callback to the settings */
    private _connectSettings() {
        // ======= x axis settings
        this._settings.x.select.onchange = () => {
            this._current.x = this._settings.x.select.value;
            const values = this._xValues();
            this._restyle({ x: values }, [0, 1]);
            this._relayout({
                'xaxis.title': this._current.x,
                'scene.xaxis.title': this._current.x,
            } as unknown as Layout);
        }

        this._settings.x.scale.onchange = () => {
            if (this._is3D()) {
                this._relayout({
                    'scene.xaxis.type': this._settings.x.scale.value,
                } as unknown as Layout);
            } else {
                this._relayout({ 'xaxis.type': this._settings.x.scale.value as Plotly.AxisType });
            }
        }

        const xRangeChange = () => {
            const min = parseFloat(this._settings.x.min.value);
            const max = parseFloat(this._settings.x.max.value);
            if (this._is3D()) {
                this._relayout({
                    'scene.xaxis.range': [min, max],
                } as unknown as Layout);
            } else {
                this._relayout({ 'xaxis.range': [min, max] });
            }
        }
        this._settings.x.min.onchange = xRangeChange;
        this._settings.x.max.onchange = xRangeChange;

        // ======= y axis settings
        this._settings.y.select.onchange = () => {
            this._current.y = this._settings.y.select.value;
            const values = this._yValues();
            this._restyle({ y:  values}, [0, 1]);
            this._relayout({
                'yaxis.title': this._current.y,
                'scene.yaxis.title': this._current.y,
            } as unknown as Layout);
        }

        this._settings.y.scale.onchange = () => {
            if (this._is3D()) {
                this._relayout({
                    'scene.yaxis.type': this._settings.y.scale.value,
                } as unknown as Layout);
            } else {
                this._relayout({ 'yaxis.type': this._settings.y.scale.value as Plotly.AxisType });
            }
        }

        const yRangeChange = () => {
            const min = parseFloat(this._settings.y.min.value);
            const max = parseFloat(this._settings.y.max.value);
            if (this._is3D()) {
                this._relayout({
                    'scene.yaxis.range': [min, max],
                } as unknown as Layout);
            } else {
                this._relayout({ 'yaxis.range': [min, max] });
            }
        }
        this._settings.y.min.onchange = yRangeChange;
        this._settings.y.max.onchange = yRangeChange;

        // ======= z axis settings
        this._settings.z.select.onchange = () => {
            if (this._settings.z.select.value === "") {
                const was3D = this._is3D();
                this._current.z = undefined;
                if (was3D) {
                    this._switch2D();
                }
            } else {
                const was3D = this._is3D();
                this._current.z = this._settings.z.select.value;
                if (!was3D) {
                    this._switch3D();
                }
            }

            const values = this._zValues();
            this._restyle({ z: values } as Data, [0, 1]);
            this._relayout({
                'scene.zaxis.title': this._current.z,
            } as unknown as Layout);
        }

        this._settings.z.scale.onchange = () => {
            this._relayout({
                'scene.zaxis.type': this._settings.z.scale.value,
            } as unknown as Layout);
        }

        const zRangeChange = () => {
            const min = parseFloat(this._settings.z.min.value);
            const max = parseFloat(this._settings.z.max.value);
            this._relayout({
                'scene.zaxis.range': [min, max],
            } as unknown as Layout);
        }
        this._settings.z.min.onchange = zRangeChange;
        this._settings.z.max.onchange = zRangeChange;

        // ======= color axis settings
        this._settings.color.select.onchange = () => {
            if (this._settings.color.select.value !== "") {
                this._current.color = this._settings.color.select.value;
                this._settings.color.min.disabled = false;
                this._settings.color.max.disabled = false;
                this._settings.colorReset.disabled = false;

                const values = this._colors(0)[0] as number[];
                const min = Math.min(...values);
                const max = Math.max(...values);

                this._settings.color.min.value = min.toString();
                this._settings.color.max.value = max.toString();

                this._relayout({
                    'coloraxis.colorbar.title.text': this._current.color,
                    'coloraxis.showscale': true,
                    'coloraxis.cmin': min,
                    'coloraxis.cmax': max,
                    'coloraxis2.cmin': min,
                    'coloraxis2.cmax': max,
                } as unknown as Layout);

            } else {
                this._current.color = undefined;
                this._settings.color.min.disabled = true;
                this._settings.color.max.disabled = true;
                this._settings.colorReset.disabled = true;
                this._settings.color.min.value = "0";
                this._settings.color.max.value = "0";

                this._relayout({
                    'coloraxis.colorbar.title.text': undefined,
                    'coloraxis.showscale': false,
                } as unknown as Layout);
            }

            this._restyle({
                'hovertemplate': this._hovertemplate(),
                'marker.color': this._colors(0),
                'marker.line.color': this._lineColors(0),
            } as Data, 0);
        }

        const colorRangeChange = () => {
            const min = parseFloat(this._settings.color.min.value);
            const max = parseFloat(this._settings.color.max.value);
            this._relayout({
                'coloraxis.cmin': min,
                'coloraxis.cmax': max,
                'coloraxis2.cmin': min,
                'coloraxis2.cmax': max,
            } as unknown as Layout);
        }
        this._settings.color.min.onchange = colorRangeChange;
        this._settings.color.max.onchange = colorRangeChange;

        this._settings.colorReset.onclick = () => {
            const values = this._colors(0)[0] as number[];
            const min = Math.min(...values);
            const max = Math.max(...values);
            this._settings.color.min.value = min.toString();
            this._settings.color.max.value = max.toString();
            this._relayout({
                'coloraxis.cmin': min,
                'coloraxis.cmax': max,
                'coloraxis2.cmin': min,
                'coloraxis2.cmax': max,
            } as unknown as Layout);
        }

        // ======= color palette
        this._settings.palette.onchange = () => {
            this._current.colorscale = this._settings.palette.value;
            this._relayout({
                'coloraxis.colorscale': this._colorScale(),
                'coloraxis2.colorscale': this._lineColorScale(),
            } as unknown as Layout);
        }

        // ======= markers symbols
        this._settings.symbol.onchange = () => {
            if (this._settings.symbol.value !== "") {
                this._current.symbols = this._settings.symbol.value;
            } else {
                this._current.symbols = undefined;
            }

            this._restyle({ 'marker.symbol': this._symbols() }, [0, 1]);

            this._restyle({
                'showlegend': this._showlegend(),
                'name': this._legendNames(),
            } as unknown as Data);

            this._relayout({
                'coloraxis.colorbar.len': this._colorbarLen(),
            } as unknown as Layout);
        }

        // ======= markers size
        this._settings.size.select.onchange = () => {
            if (this._settings.size.select.value != "") {
                this._current.size = this._settings.size.select.value;
            } else {
                this._current.size = undefined;
            }

            const factor = parseInt(this._settings.size.factor.value);
            this._restyle({ 'marker.size': this._sizes(factor, 0) } as Data, 0);
        };

        this._settings.size.factor.onchange = () => {
            const factor = parseInt(this._settings.size.factor.value);
            this._restyle({ 'marker.size': this._sizes(factor, 0) } as Data, 0);
        }
    }

    /**
     * Fill possible values for the settings, depending on the properties in
     * the dataset
     */
    private _setupSettings() {
        // ============== Setup the map options ==============
        // ======= data used as x values
        this._settings.x.select.options.length = 0;
        for (const key in this._properties()) {
            this._settings.x.select.options.add(new Option(key, key));
        }
        this._settings.x.select.selectedIndex = 0;

        // ======= data used as y values
        this._settings.y.select.options.length = 0;
        for (const key in this._properties()) {
            this._settings.y.select.options.add(new Option(key, key));
        }
        this._settings.y.select.selectedIndex = 1;

        // ======= data used as z values
        // first option is "none"
        this._settings.z.select.options.length = 1;
        for (const key in this._properties()) {
            this._settings.z.select.options.add(new Option(key, key));
        }
        this._settings.z.select.selectedIndex = 0;

        // ======= data used as color values
        this._settings.color.select.options.length = 1;
        for (const key in this._properties()) {
            this._settings.color.select.options.add(new Option(key, key));
        }
        if (this._hasColors()) {
            // index 0 is 'none', 1 is the x values, 2 the y values, use 3 for
            // colors
            this._settings.color.select.selectedIndex = 3;
            this._settings.color.min.disabled = false;
            this._settings.color.max.disabled = false;
            this._settings.colorReset.disabled = false;

            const values = this._colors(0)[0] as number[];
            const min = Math.min(...values);
            const max = Math.max(...values);

            this._settings.color.min.value = min.toString();
            this._settings.color.max.value = max.toString();
        } else {
            this._settings.color.select.selectedIndex = 0;
            this._settings.color.min.disabled = true;
            this._settings.color.max.disabled = true;
            this._settings.colorReset.disabled = true;

            this._settings.color.min.value = "0";
            this._settings.color.max.value = "0";
        }

        // ======= color palette
        this._settings.palette.options.length = 0;
        for (const key in COLOR_MAPS) {
            this._settings.palette.options.add(new Option(key, key));
        }
        this._settings.palette.value = this._current.colorscale;

        // ======= marker symbols
        // first option is "default"
        this._settings.symbol.options.length = 1;
        this._settings.symbol.selectedIndex = 0;
        for (const key in this._properties()) {
            if (this._properties()[key].string !== undefined) {
                this._settings.symbol.options.add(new Option(key, key));
            }
        }

        // ======= marker size
        // first option is "default"
        this._settings.size.select.options.length = 1;
        this._settings.size.select.selectedIndex = 0;
        for (const key in this._properties()) {
            this._settings.size.select.options.add(new Option(key, key));
        }
        this._settings.size.factor.value = "75";
    }

    /** Actually create the Plotly plot */
    private _createPlot() {
        this._plot.innerHTML = '';

        const colors = this._colors();
        const lineColors = this._lineColors();
        // default value for the size factor is 75
        const sizes = this._sizes(75);
        const symbols = this._symbols();

        const x = this._xValues();
        const y = this._yValues();
        const z = this._zValues();

        // The main trace, containing default data
        const main = {
            type: "scattergl",
            name: this._name,
            x: x[0],
            y: y[0],
            z: z[0],
            hovertemplate: this._hovertemplate(),
            showlegend: false,
            mode: "markers",
            marker: {
                color: colors[0],
                size: sizes[0],
                symbol: symbols[0],
                line: {
                    color: lineColors[0],
                    width: 1,
                    coloraxis: "coloraxis2",
                },
                coloraxis: "coloraxis",
                // prevent plolty from messing with opacity when doing bubble
                // style charts (different sizes for each point)
                opacity: 1,
            },
        };

        // Create a second trace to store the last clicked point, in order to
        // display it on top of the main plot with different styling. This is
        // only used in 3D mode, since it is way slower than moving
        // this._selectedMarker around.
        const selected = {
            type: "scattergl",
            name: "",
            x: [NaN],
            y: [NaN],
            z: [NaN],
            hoverinfo: "none",
            showlegend: false,
            mode: "markers",
            marker: {
                color: colors[1],
                size: sizes[1],
                line: {
                    color: lineColors[1],
                    width: 0.5,
                },
            },
        };
        const traces = [main as Data, selected as Data];

        // add empty traces to be able to display the symbols legend
        // one trace for each possible symbol
        for (let i=0; i<this._data.maxSymbols; i++) {
            const data = {
                type: "scattergl",
                name: "",
                x: [NaN],
                y: [NaN],
                z: [NaN],
                showlegend: false,
                mode: "markers",
                marker: {
                    color: "black",
                    size: sizes[0],
                    symbol: i,
                },
            };
            traces.push(data as Data);
        }

        let layout = {
            ...DEFAULT_LAYOUT
        };
        layout.title.text = this._name;
        layout.xaxis.title = this._current.x;
        layout.yaxis.title = this._current.y;
        layout.xaxis.type = this._settings.x.scale.value;
        layout.yaxis.type = this._settings.y.scale.value;
        layout.scene.xaxis.title = this._current.x;
        layout.scene.yaxis.title = this._current.y;
        layout.scene.zaxis.title = this._current.z;
        layout.coloraxis.colorscale = this._colorScale();
        layout.coloraxis.cmin = parseFloat(this._settings.color.min.value);
        layout.coloraxis.cmax = parseFloat(this._settings.color.max.value);
        layout.coloraxis2.cmin = parseFloat(this._settings.color.min.value);
        layout.coloraxis2.cmax = parseFloat(this._settings.color.max.value);
        layout.coloraxis2.colorscale = this._lineColorScale();
        layout.coloraxis.colorbar.title.text = this._current.color;

        // Create an empty plot and fill it below
        Plotly.newPlot(this._plot, traces, layout as Partial<Layout>, DEFAULT_CONFIG as Config)
            .catch(e => setTimeout(() => {throw e;}));

        this._plot.on("plotly_click", (event: Plotly.PlotMouseEvent) => {
            const environment = event.points[0].pointNumber;
            const indexes = this._indexer.from_environment(environment);
            this.select(indexes)
            this.onselect(indexes);
        });
        this._plot.on("plotly_afterplot", () => this._afterplot());

        this._updateSelectedMarker();
    }

    /** Get the currently available properties: either `"atom"` or `"structure"` properties */
    private _properties(): {[name: string]: NumericProperty} {
        return this._data[this._indexer.mode]
    }

    /** Get the plotly hovertemplate depending on `this._current.color` */
    private _hovertemplate(): string {
        if (this._hasColors()) {
            return this._current.color + ": %{marker.color:.2f}<extra></extra>";
        } else {
            return "%{x:.2f}, %{y:.2f}<extra></extra>";
        }
    }

    /**
     * Get the values to use for the x axis with the given plotly `trace`,
     * or all of them if `trace === undefined`
     */
    private _xValues(trace?: number): Array<number[]> {
        const values = this._properties()[this._current.x].values;
        return this._selectTrace(values, [values[this._selected]], trace);
    }

    /**
     * Get the values to use for the y axis with the given plotly `trace`,
     * or all of them if `trace === undefined`
     */
    private _yValues(trace?: number): Array<number[]> {
        const values = this._properties()[this._current.y].values;
        return this._selectTrace(values, [values[this._selected]], trace);
    }

    /**
     * Get the values to use for the z axis with the given plotly `trace`,
     * or all of them if `trace === undefined`
     */
    private _zValues(trace?: number): Array<undefined | number[]> {
        if (!this._is3D()) {
            return this._selectTrace(undefined, undefined, trace);
        }

        const values = this._properties()[this._current.z!].values;
        return this._selectTrace(values, [values[this._selected]], trace);
    }

    /**
     * Get the color values to use with the given plotly `trace`, or all of
     * them if `trace === undefined`
     */
    private _colors(trace?: number): Array<string | number | number[]> {
        let values;
        if (this._hasColors()) {
            values = this._properties()[this._current.color!].values;
        } else {
            values = 0.5;
        }

        return this._selectTrace<string | number | number[]>(values, "#007bff", trace);
    }

    /**
     * Get the **line** color values to use with the given plotly `trace`, or
     * all of them if `trace === undefined`
     */
    private _lineColors(trace?: number): Array<string | number | number[]> {
        let values;
        if (this._hasColors()) {
            values = this._properties()[this._current.color!].values;
        } else {
            values = 0.5;
        }

        return this._selectTrace<string | number | number[]>(values, "black", trace);
    }

    /** Get the colorscale to use for markers in the main plotly trace */
    private _colorScale(): Plotly.ColorScale {
        if (this._is3D()) {
            return COLOR_MAPS[this._current.colorscale].rgb;
        } else {
            return COLOR_MAPS[this._current.colorscale].rgba;
        }
    }

    /** Get the colorscale to use for markers lines in the main plotly trace */
    private _lineColorScale(): Plotly.ColorScale {
        return COLOR_MAPS[this._current.colorscale].rgb;
    }

    /**
     * Get the values to use as marker size with the given plotly `trace`, or
     * all of them if `trace === undefined`.
     *
     * The size scaling parameter should be given in `sizeSliderValue`.
     */
    private _sizes(sizeSliderValue: number, trace?: number): Array<number | number[]> {
        // Transform the linear value from the slider into a logarithmic scale
        const logSlider = (value: number) => {
            const min_slider = 1;
            const max_slider = 100;

            const min_value = Math.log(1.0 / 6.0);
            const max_value = Math.log(2.0);

            const factor = (max_value - min_value) / (max_slider - min_slider);
            return Math.exp(min_value + factor * (value - min_slider));
        }

        const factor = logSlider(sizeSliderValue);

        const defaultSize = this._is3D() ? 4.5 : 10;
        let values;
        if (this._current.size === undefined) {
            values = defaultSize * factor;
        } else {
            const sizes = this._properties()[this._current.size].values;
            const min = Math.min.apply(Math, sizes);
            const max = Math.max.apply(Math, sizes);
            const defaultSize = this._is3D() ? 12 : 20;
            // normalize inside [0, 10 * factor]
            values = sizes.map((v: number) => {
                const scaled = (v - min) / (max - min);
                return defaultSize * factor * (scaled + 0.05);
            })
        }

        return this._selectTrace(values, 2 * defaultSize, trace);
    }

    /**
     * Get the values to use as marker symbol with the given plotly `trace`, or
     * all of them if `trace === undefined`.
     */
    private _symbols(trace?: number): Array<number | number[]> {
        if (this._current.symbols !== undefined) {
            assert(!this._is3D());
            const values = this._properties()[this._current.symbols].values;
            return this._selectTrace<number | number[]>(values, values[this._selected], trace);
        } else {
            // default to 0 (i.e. circles)
            return this._selectTrace<number | number[]>(0, 0, trace);
        }
    }

    /** Should we show the legend for the various symbols used? */
    private _showlegend(): Array<boolean> {
        const result = [false, false];

        if (this._current.symbols !== undefined) {
            assert(!this._is3D());
            for (let i=0; i<this._symbolsCount()!; i++) {
                result.push(true);
            }
            return result;
        } else {
            for (let i=0; i<this._data.maxSymbols; i++) {
                result.push(false);
            }
            return result;
        }
    }

    /** Get the list of symbols names to use for the legend */
    private _legendNames(): Array<string> {
        const result = [this._name, ""];

        if (this._current.symbols !== undefined) {
            assert(!this._is3D());
            const names = this._properties()[this._current.symbols].string!.strings();
            for (const name of names) {
                result.push(name);
            }
            return result;
        } else {
            for (let i=0; i<this._data.maxSymbols; i++) {
                result.push("");
            }
            return result;
        }
    }

    /**
     * Select either main, selected or both depending on `trace`, and return
     * them in a mode usable with `Plotly.restyle`/[[PropertiesMap._restyle]]
     */
    private _selectTrace<T>(main: T, selected: T, trace?: number): Array<T> {
        if (trace === 0) {
            return [main];
        } else if (trace === 1) {
            return [selected]
        } else if (trace === undefined) {
            return [main, selected]
        } else {
            throw Error("internal error: invalid trace number")
        }
    }

    /** How many different symbols are being displayed */
    private _symbolsCount(): number {
        if (this._current.symbols !== undefined) {
            return this._properties()[this._current.symbols].string!.strings().length;
        } else {
            return 0;
        }
    }

    /** Get the length of the colorbar to accomodate for the legend */
    private _colorbarLen(): number {
        /// Heigh of a legend item in plot unit
        const LEGEND_ITEM_HEIGH = 0.045;
        return 1 - LEGEND_ITEM_HEIGH * this._symbolsCount()
    }

    /** Does the current plot use color values? */
    private _hasColors(): boolean {
        return this._current.color !== undefined;
    }

    /** Is the the current plot a 3D plot? */
    private _is3D(): boolean {
        return this._current.z !== undefined;
    }

    /** Switch current plot from 2D to 3D */
    private _switch3D() {
        assert(this._is3D());
        // Can not use symbols with 3D plots, why?
        this._settings.symbol.disabled = true;
        this._current.symbols = undefined;
        this._settings.z.scale.disabled = false;
        this._settings.z.min.disabled = false;
        this._settings.z.max.disabled = false;

        // Switch all traces to 3D mode
        this._restyle({
            type: 'scatter3d',
            showlegend: this._showlegend(),
            name: this._legendNames(),
        } as unknown as Data);

        this._selectedMarker.style.display = "none";
        this._updateSelectedMarker();

        // Change the data that vary between 2D and 3D mode
        const factor = parseInt(this._settings.size.factor.value);
        this._restyle({
            // size and symbols change from 2D to 3D
            'marker.size': this._sizes(factor),
            'marker.symbol': this._symbols(),
            'marker.line.width': [0, 1],
        } as Data, [0, 1]);

        this._relayout({
            // Do not use opacity in 3D mode, since it renders horribly
            'coloraxis.colorscale': this._colorScale(),
            // change colorbar length to accomodate for symbols legend
            'coloraxis.colorbar.len': this._colorbarLen(),
            // Carry over axis types
            'scene.xaxis.type': this._settings.x.scale.value as Plotly.AxisType,
            'scene.yaxis.type': this._settings.y.scale.value as Plotly.AxisType,
            'scene.zaxis.type': this._settings.z.scale.value as Plotly.AxisType,
        } as unknown as Layout);
    }

    /** Switch current plot from 3D back to 2D */
    private _switch2D() {
        assert(!this._is3D());
        this._settings.symbol.disabled = false;
        if (this._settings.symbol.value !== "") {
            this._current.symbols = this._settings.symbol.value;
        }
        this._settings.z.scale.disabled = true;
        this._settings.z.min.disabled = true;
        this._settings.z.max.disabled = true;

        // switch all traces to 2D mode
        this._restyle({
            type: 'scattergl',
            showlegend: this._showlegend(),
            name: this._legendNames(),
        } as unknown as Data);

        this._selectedMarker.style.display = "block";
        this._restyle({x: [NaN], y: [NaN]}, 1);
        this._updateSelectedMarker();

        // Change the data that vary between 2D and 3D mode
        const factor = parseInt(this._settings.size.factor.value);
        this._restyle({
            // size and symbols change from 2D to 3D
            'marker.size': this._sizes(factor),
            'marker.symbol': this._symbols(),
            'marker.line.width': [1, 1],
        } as Data, [0, 1]);

        this._relayout({
            // Do not use opacity in 3D mode, since it renders horribly
            'coloraxis.colorscale': this._colorScale(),
            // change colorbar length to accomodate for symbols legend
            'coloraxis.colorbar.len': this._colorbarLen(),
            // Carry over axis types
            'xaxis.type': this._settings.x.scale.value as Plotly.AxisType,
            'yaxis.type': this._settings.y.scale.value as Plotly.AxisType,
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
            this._settings.x.min.value = layout.scene.xaxis.range[0].toString();
            this._settings.x.max.value = layout.scene.xaxis.range[1].toString();

            this._settings.y.min.value = layout.scene.yaxis.range[0].toString();
            this._settings.y.max.value = layout.scene.yaxis.range[1].toString();

            this._settings.z.min.value = layout.scene.zaxis.range[0].toString();
            this._settings.z.max.value = layout.scene.zaxis.range[1].toString();
        } else {
            this._settings.x.min.value = layout.xaxis.range[0].toString();
            this._settings.x.max.value = layout.xaxis.range[1].toString();

            this._settings.y.min.value = layout.yaxis.range[0].toString();
            this._settings.y.max.value = layout.yaxis.range[1].toString();

            this._updateSelectedMarker();
        }
    }

    /** Update the position of the selected marker */
    private _updateSelectedMarker() {
        if (this._is3D()) {
            this._restyle({
                x: this._xValues(1),
                y: this._yValues(1),
                z: this._zValues(1),
            } as Data, 1);
        } else {
            const xaxis = this._plot._fullLayout.xaxis;
            const yaxis = this._plot._fullLayout.yaxis;

            const computeX = (data: number) => xaxis.l2p(data) + xaxis._offset;
            const computeY = (data: number) => yaxis.l2p(data) + yaxis._offset;

            const x = computeX(this._xValues(0)[0][this._selected]);
            const y = computeY(this._yValues(0)[0][this._selected]);

            // hide the point if it is outside the plot, allow for up to 10px
            // overflow (for points just on the border)
            const minX = computeX(xaxis.range[0]) - 10;
            const maxX = computeX(xaxis.range[1]) + 10;
            const minY = computeY(yaxis.range[1]) - 10;
            const maxY = computeY(yaxis.range[0]) + 10;
            if (x < minX || x > maxX || y < minY || y > maxY) {
                this._selectedMarker.style.display = "none";
            } else {
                this._selectedMarker.style.display = "block";
            }

            // -10 since we want the centers to match, and the marker div is 20px wide
            this._selectedMarker.style.top = `${y - 10}px`;
            this._selectedMarker.style.left = `${x - 10}px`;
        }
    }
}
