import * as Plotly from "./lib/plotly-scatter";
import {Config, Data, Layout, PlotlyHTMLElement} from "./lib/plotly-scatter";

import {COLOR_MAPS} from "./colorscales";
import {make_draggable} from "./draggable";

import {MapInput, MapData} from "./map_data"

const HTML_SETTINGS = require("./static/settings.html");

const DEFAULT_LAYOUT = {
    hovermode: "closest",
    showlegend: false,
    xaxis: {
        title: "",
        zeroline: false,
    },
    yaxis: {
        title: "",
        zeroline: false,
    },
    margin: {
        l: 50,
        t: 50,
        r: 50,
        b: 50,
    }
};

const DEFAULT_CONFIG = {
    displayModeBar: true,
    displaylogo: false,
    scrollZoom: true,
    modeBarButtonsToRemove: [
        "hoverClosestCartesian",
        "hoverCompareCartesian",
        "toggleSpikelines",
        "resetScale2d",
        "select2d",
        "lasso2d",
    ],
};

function getByID<HTMLType>(id: string): HTMLType {
    const e = document.getElementById(id);
    if (e === null) {
        throw Error(`unable to get element with id ${id}`);
    }
    return e as unknown as HTMLType;
}

export class ScatterPlot {
    /// HTML root holding the full plot
    private _root: HTMLElement;
    /// Plotly plot
    private _plot!: PlotlyHTMLElement;
    /// The dataset name
    private _name: string;
    /// All known properties
    private _data: MapData
    /// Storing the callback for when the plot is clicked
    private _selectedCallback: (index: number) => void;
    /// Index of the currently selected point
    private _selected: number;
    /// Currently displayed data
    private _current!: {
        /// Name of the properties in `this._data` used for x values
        x: string,
        /// Name of the properties in `this._data` used for y values
        y: string,
        /// Name of the properties in `this._data` used for color values
        color: string,
        /// Symbol used for all marker / every marker
        symbols: number | number[],
    }
    /// callback to get the initial positioning of the settings modal. The
    /// callback gets the current placement of the settings as a DOMRect, and
    /// should return top and left positions in pixels, used with
    /// `position: fixed`
    private _settingsPlacement!: (rect: DOMRect) => {top: number, left: number};

    constructor(id: string, name: string, properties: MapInput) {
        this._name = name;
        this._selectedCallback = (_) => { return; };
        this._selected = 0;

        const root = document.getElementById(id);
        if (root === null) {
            throw Error(`could not find HTML element #${id}`)
        }
        this._root = root;
        this._root.style.position = 'relative';

        this._data = new MapData(properties);
        this._setupDefaults();

        this._createSettings();
        this._setupSettings();

        this._createPlot();
        this._setupPlot();
    }

    /// Register a callback to be called when the selected envirronement is
    /// updated
    public onSelectedUpdate(callback: (index: number) => void) {
        this._selectedCallback = callback;
    }

    /// Change the selected environement to the one with the given `index`
    public select(index: number) {
        this._selected = index;

        let symbol;
        if (typeof(this._current.symbols) === "number") {
            symbol = this._current.symbols;
        } else {
            symbol = this._current.symbols[this._selected];
        }

        Plotly.restyle(this._plot, {
            x: [[this._data[this._current.x].values[this._selected]]],
            y: [[this._data[this._current.y].values[this._selected]]],
            "marker.symbol": [symbol],
        }, 1);
        this._selectedCallback(this._selected);
    }

    /// Change the current dataset to the provided one, without re-creating the
    /// plot
    public changeDataset(name: string, properties: MapInput) {
        this._name = name;
        this._selected = 0;
        this._data = new MapData(properties);
        this._setupDefaults();
        this._setupSettings();
        this._setupPlot();
    }

    /// Use the given callback to compute the placement of the settings modal.
    ///
    /// The callback gets the current placement of the settings as a DOMRect,
    /// and should return top and left positions in pixels, used with `position:
    /// fixed`. The callback is called once, the first time the settings are
    /// opened.
    public computeSettingsPlacement(callback: (rect: DOMRect) => {top: number, left: number}) {
        this._settingsPlacement = callback;
    }

    private _setupDefaults() {
        const prop_names = Object.keys(this._data);
        if (prop_names.length < 2) {
            throw Error("Sketchmap needs at least two properties to plot")
        }
        this._current = {
            x: prop_names[0],
            y: prop_names[1],
            color: (prop_names.length > 2) ? prop_names[2] : "",
            symbols: 0
        }
    }

    private _createSettings() {
        // use HTML5 template to generate a DOM object from an HTML string
        const template = document.createElement('template');
        template.innerHTML = `<button
            href="#skv-settings"
            class="btn btn-light btn-sm skv-open-settings"
            data-toggle="modal">
                <div class="skv-hamburger"><div></div><div></div><div></div></div>
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
        make_draggable(modalDialog, ".modal-header");

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

    private _setupSettings() {
        // ============== Setup the map options ==============
        // ======= data used as x values
        const xValues = getByID<HTMLSelectElement>('skv-x');
        xValues.options.length = 0;
        for (const key in this._data) {
            xValues.options.add(new Option(key, key));
        }
        xValues.selectedIndex = 0;
        xValues.onchange = () => {
            this._current.x = xValues.value;
            const values = this._data[xValues.value].values;
            const data = {
                x: [values, [values[this._selected]]]
            }
            Plotly.restyle(this._plot, data).catch(e => console.error(e));
            Plotly.relayout(this._plot, {
                'xaxis.title': xValues.value
            }).catch(e => console.error(e));
        }

        // ======= data used as y values
        const yValues = getByID<HTMLSelectElement>('skv-y');
        yValues.options.length = 0;
        for (const key in this._data) {
            yValues.options.add(new Option(key, key));
        }
        yValues.selectedIndex = 1;
        yValues.onchange = () => {
            this._current.y = yValues.value;
            const values = this._data[yValues.value].values;
            const data = {
                y: [values, [values[this._selected]]]
            }
            Plotly.restyle(this._plot, data).catch(e => console.error(e));
            Plotly.relayout(this._plot, {
                'yaxis.title': yValues.value
            }).catch(e => console.error(e));
        }

        // ======= marker color
        const color = getByID<HTMLSelectElement>('skv-color');
        color.options.length = 1;
        for (const key in this._data) {
            color.options.add(new Option(key, key));
        }
        if (this._current.color !== "") {
            // 2 + 1 since the first item is 'none'
            color.selectedIndex = 3;
        }
        color.onchange = () => {
            this._current.color = color.value;
            const values = (this._current.color === "") ? [0.5] : [this._data[this._current.color].values];
            const data = {
                hovertemplate: this._hovertemplate(),
                'marker.color': values,
                'marker.line.color': values,
                'marker.colorbar.title': this._current.color,
            };
            Plotly.restyle(this._plot, data, 0).catch(e => console.error(e));
        }

        // ======= color palette
        const palette = getByID<HTMLSelectElement>('skv-palette');
        palette.options.length = 0;
        for (const key in COLOR_MAPS) {
            palette.options.add(new Option(key, key));
        }
        palette.value = 'inferno';

        palette.onchange = () => {
            const data = {
                'marker.colorscale': [COLOR_MAPS[palette.value].rgba],
                'marker.line.colorscale': [COLOR_MAPS[palette.value].rgb],
            };
            Plotly.restyle(this._plot, data, 0).catch(e => console.error(e));
        }

        // ======= marker symbols
        const symbol = getByID<HTMLSelectElement>('skv-symbol');
        symbol.options.length = 1;
        for (const key in this._data) {
            if (this._data[key].string !== undefined) {
                symbol.options.add(new Option(key, key));
            }
        }

        symbol.onchange = () => {
            let selected;
            if (symbol.value !== "") {
                this._current.symbols = this._data[symbol.value].values;
                selected = this._current.symbols[this._selected];
            } else {
                // reset default
                this._current.symbols = 0;
                selected = 0;
            }
            const data = {
                'marker.symbol': [this._current.symbols, selected],
            };
            Plotly.restyle(this._plot, data).catch(e => console.error(e));
        }

        // ======= marker size
        const size = getByID<HTMLSelectElement>('skv-size');
        size.options.length = 1;
        const sizeFactor = getByID<HTMLInputElement>('skv-size-factor');
        sizeFactor.value = "25";
        for (const key in this._data) {
            size.options.add(new Option(key, key));
        }

        const changeSize = () => {
            // Transform the linear value from the slider into a logarithmic scale
            const logSlider = (value: number) => {
                const min_slider = parseInt(sizeFactor.min);
                const max_slider = parseInt(sizeFactor.max);

                const min_value = Math.log(1.0 / 6.0);
                const max_value = Math.log(6.0);

                const factor = (max_value - min_value) / (max_slider - min_slider);
                return Math.exp(min_value + factor * (value - min_slider));
            }

            const factor = logSlider(parseInt(sizeFactor.value));
            let markerSize;

            if (size.value === "") {
                markerSize = [10 * factor];
            } else {
                const values = this._data[size.value].values;
                const min = Math.min.apply(Math, values);
                const max = Math.max.apply(Math, values);
                markerSize = [values.map((v) => {
                    const scaled = (v - min) / (max - min);
                    return 10 * factor * (scaled + 0.05);
                })]
            }

            const data = {
                'marker.size': markerSize,
            };
            Plotly.restyle(this._plot, data, 0).catch(e => console.error(e));
        }
        size.onchange = changeSize;
        sizeFactor.onchange = changeSize;
    }

    private _createPlot() {
        this._plot = document.createElement("div") as unknown as PlotlyHTMLElement;
        this._plot.style.width = "100%";
        this._plot.style.height = "100%";
        this._plot.style.minHeight = "550px";
        this._root.appendChild(this._plot);

        // create an empty plot, and fill it in _setupPlot
        Plotly.newPlot(this._plot, [{}, {}], {}, DEFAULT_CONFIG as Config)
            .catch(e => console.error(e));

        this._plot.on("plotly_click", (event: Plotly.PlotMouseEvent) => this.select(event.points[0].pointNumber));
    }

    private _setupPlot() {
        const color = (this._current.color === "") ? 0.5 : this._data[this._current.color].values;
        const fullData = {
            type: "scattergl",
            x: this._data[this._current.x].values,
            y: this._data[this._current.y].values,
            hovertemplate: this._hovertemplate(),
            mode: "markers",
            marker: {
                color: color,
                colorscale: COLOR_MAPS.inferno.rgba,
                line: {
                    color: color,
                    colorscale: COLOR_MAPS.inferno.rgb,
                    width: 1,
                },
                size: 10,
                symbol: this._current.symbols,
                showscale: true,
                colorbar: {
                    title: this._current.color,
                    thickness: 20,
                }
            },
        };

        // Create a second trace to store the last clicked point, in order to
        // display it on top of the main plot with different styling
        const selected = {
            type: "scattergl",
            x: [fullData.x[this._selected]],
            y: [fullData.y[this._selected]],
            mode: "markers",
            marker: {
                color: "#007bff",
                line: {
                    color: "black",
                    width: 0.5,
                },
                size: 20,
                symbol: this._current.symbols,
            },
            hoverinfo: "none",
        };

        const layout = {
            ...DEFAULT_LAYOUT,
            title: {
                text: this._name,
                font: {
                    size: 25,
                },
            },
        };

        layout.xaxis.title = this._current.x;
        layout.yaxis.title = this._current.y;

        Plotly.react(
            this._plot, [fullData as Data, selected as Data],
            layout as Layout,
        ).catch(e => console.error(e));
        this._plot.on("plotly_click", (event: Plotly.PlotMouseEvent) => this.select(event.points[0].pointNumber));
    }

    /// Get the plotly hovertemplate depending on `this._current.color`
    private _hovertemplate(): string {
        if (this._current.color !== "") {
            return this._current.color + ": %{marker.color:.2f}<extra></extra>";
        } else {
            return "%{x:.2f}, %{y:.2f}<extra></extra>";
        }
    }
}
