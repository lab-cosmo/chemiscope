import * as Plotly from "plotly.js";
import {Config, Data, Layout, PlotlyHTMLElement} from "plotly.js";
import {COLOR_MAPS} from "./colorscales";
import {make_draggable} from "./draggable";

require('./static/sketchviz.css');
const HTML_SETTINGS = require("./static/settings.html");

interface MapData {
    [key: string]: number[] | string[];
}

interface MapInput {
    name: string;
    data: MapData;
}

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

export class Sketchmap {
    /// HTML root holding the full plot
    private _root: HTMLElement;
    /// Plotly plot
    private _plot!: PlotlyHTMLElement;
    /// The dataset name
    private _name: string;
    /// All known properties
    private _data: {
        numeric: {
            [key: string]: number[]
        };
        strings: {
            [key: string]: string[]
        };
    };
    /// Storing the callback for when the plot is clicked
    private _selectedCallback: (index: number) => void;
    /// Index of the currently selected point
    private _selected: number;
    /// Currently displayed data, in this._data
    private _current!: {
        x: string,
        y: string,
        color?: string,
    }

    constructor(id: string, data: MapInput) {
        this._name = data.name;
        this._selectedCallback = (_) => { return; };
        this._selected = 0;
        this._data = {
            numeric: {},
            strings: {},
        };

        const root = document.getElementById(id);
        if (root === null) {
            throw Error(`could not find HTML element #${id}`)
        }
        this._root = root;
        this._root.classList.add('skv-root');

        this._extractProperties(data.data);
        this._setupSettings();
        this._createPlot();
    }

    /// Register a callback to be called when the selected envirronement is
    /// updated
    public onSelectedUpdate(callback: (index: number) => void) {
        this._selectedCallback = callback;
    }

    /// Change the selected environement to the one with the given `index`
    public select(index: number) {
        this._selected = index;
        Plotly.restyle(this._plot, {
            x: [[this._data.numeric[this._current.x][this._selected]]],
            y: [[this._data.numeric[this._current.y][this._selected]]],
        }, 1);
        this._selectedCallback(this._selected);
    }

    private _extractProperties(data: MapData) {
        // check that all properties have the same size
        let size = -1;
        for (const key in data) {
            if (size === -1) {
                size = data[key].length;
            }

            if (data[key].length !== size) {
                throw Error("not all properties have the same size")
            }
        }

        if (size === 0) {
            return;
        }

        for (const key in data) {
            const prop_type = typeof(data[key][0]);
            if (prop_type === "number") {
                this._data.numeric[key] = data[key] as number[];
            } else if (prop_type === "string") {
                this._data.strings[key] = data[key] as string[];
            } else {
                throw Error(`unexpected property type '${prop_type}'`);
            }
        }

        const num_prop_names = Object.keys(this._data.numeric);
        if (num_prop_names.length < 2) {
            throw Error("Sketchmap needs at least two numeric properties to plot")
        }

        this._current = {
            x: num_prop_names[0],
            y: num_prop_names[1],
        }
        if (num_prop_names.length > 2) {
            this._current.color = num_prop_names[2]
        }
    }

    private _setupSettings() {
        // use HTML5 template to generate a DOM object from an HTML string
        const template = document.createElement('template');
        template.innerHTML = `<button
            href="#skv-settings"
            class="btn btn-light btn-sm skv-open-settings"
            data-toggle="modal">
                <div class="skv-hamburger"><div></div><div></div><div></div></div>
            </button>`;
        this._root.append(template.content.firstChild!);

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

        // ============== Setup the map options ==============
        // ======= data used as x values
        const xValues = getByID<HTMLSelectElement>('skv-x');
        for (const key in this._data.numeric) {
            xValues.options.add(new Option(key, key));
        }
        xValues.selectedIndex = 0;
        xValues.onchange = () => {
            this._current.x = xValues.value;
            const data = {
                x: [
                    this._data.numeric[xValues.value],
                    [this._data.numeric[xValues.value][this._selected]]
                ]
            }
            Plotly.restyle(this._plot, data).catch(e => console.error(e));
            Plotly.relayout(this._plot, {
                'xaxis.title': xValues.value
            } as unknown as Layout).catch(e => console.error(e));
        }

        // ======= data used as y values
        const yValues = getByID<HTMLSelectElement>('skv-y');
        for (const key in this._data.numeric) {
            yValues.options.add(new Option(key, key));
        }
        yValues.selectedIndex = 1;
        yValues.onchange = () => {
            this._current.y = yValues.value;
            const data = {
                y: [
                    this._data.numeric[yValues.value],
                    [this._data.numeric[yValues.value][this._selected]]
                ]
            }
            Plotly.restyle(this._plot, data).catch(e => console.error(e));
            Plotly.relayout(this._plot, {
                'yaxis.title': yValues.value
            } as unknown as Layout).catch(e => console.error(e));
        }

        // ======= marker color
        const color = getByID<HTMLSelectElement>('skv-color');
        for (const key in this._data.numeric) {
            color.options.add(new Option(key, key));
        }
        if (this._current.color !== undefined) {
            color.selectedIndex = 2;
        }
        color.onchange = () => {
            this._current.color = color.value;
            const data = {
                'marker.color': [this._data.numeric[color.value]],
                'marker.line.color': [this._data.numeric[color.value]],
                'marker.colorbar.title': color.value,
            };
            Plotly.restyle(this._plot, data, 0).catch(e => console.error(e));
        }

        // ======= marker shapes
        const shape = getByID<HTMLSelectElement>('skv-shape');
        // TODO

        // ======= marker size
        const size = getByID<HTMLSelectElement>('skv-size');
        const sizeFactor = getByID<HTMLInputElement>('skv-size-factor');
        for (const key in this._data.numeric) {
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

            if (size.value === "default") {
                markerSize = 10 * factor;
            } else {
                const sizes = this._data.numeric[size.value];
                const min = Math.min.apply(Math, sizes);
                const max = Math.max.apply(Math, sizes);
                markerSize = sizes.map((v) => {
                    const scaled = (v - min) / (max - min);
                    return 10 * factor * (scaled + 0.05);
                })
            }

            const data = {
                'marker.size': [markerSize],
            };
            Plotly.restyle(this._plot, data as unknown as Data, 0).catch(e => console.error(e));
        }
        size.onchange = changeSize;
        sizeFactor.onchange = changeSize;

        // ======= color palette
        const palette = getByID<HTMLSelectElement>('skv-palette');
        for (const key in COLOR_MAPS) {
            palette.options.add(new Option(key, key));
        }
        palette.value = 'inferno';

        palette.onchange = () => {
            const data = {
                'marker.colorscale': [COLOR_MAPS[palette.value].rgba],
                'marker.line.colorscale': [COLOR_MAPS[palette.value].rgb],
            };
            Plotly.restyle(this._plot, data as unknown as Data, 0).catch(e => console.error(e));
        }
    }

    private _createPlot() {
        this._plot = document.createElement("div") as unknown as PlotlyHTMLElement;
        this._plot.style.width = "100%";
        this._plot.style.height = "100%";
        this._plot.style.minHeight = "550px";
        this._root.appendChild(this._plot);

        const hovertemplate = (this._current.color === undefined) ? "" :
            this._current.color + ": %{marker.color:.2f} <extra></extra>";

        const color = (this._current.color === undefined) ? 0.5 : this._data.numeric[this._current.color];
        const fullData = {
            x: this._data.numeric[this._current.x],
            y: this._data.numeric[this._current.y],
            hovertemplate: hovertemplate,
            marker: {
                color: color,
                colorscale: COLOR_MAPS.inferno.rgba,
                line: {
                    color: color,
                    colorscale: COLOR_MAPS.inferno.rgb,
                    width: 1.5,
                },
                size: 10,
                showscale: true,
                colorbar: {
                    title: this._current.color ?? "",
                    thickness: 20,
                }
            },
            mode: "markers",
            type: "scattergl",
        };

        // Create a second trace to store the last clicked point, in order to
        // display it on top of the main plot with different styling
        const clicked = {
            x: [fullData.x[this._selected]],
            y: [fullData.y[this._selected]],
            type: "scattergl",
            mode: "markers",
            marker: {
                color: "cyan",
                line: {
                    color: "black",
                    width: 0.5,
                },
                size: 18,
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

        Plotly.newPlot(
            this._plot, [fullData as Data, clicked as Data],
            layout as Layout,
            DEFAULT_CONFIG as Config,
        ).catch(e => console.error(e));
        this._plot.on("plotly_click", (event: Plotly.PlotMouseEvent) => this.select(event.points[0].pointNumber));
    }
}
