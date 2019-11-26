import * as Plotly from "plotly.js";
import {Config, Data, Layout, PlotlyHTMLElement} from "plotly.js";
import {Inferno} from "./colorscales";

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
        zeroline:false,
    },
    yaxis: {
        zeroline:false,
    },
    margin: {
        l: 50,
        t: 50,
        r: 50,
        b: 20,
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
    private _clicked_cb: (index: number) => void;
    /// Currently displayed data, in this._data
    private _current!: {
        x: string,
        y: string,
        color?: string,
    }

    constructor(id: string, data: MapInput) {
        this._name = data.name;
        this._clicked_cb = (_) => { return; };
        this._data = {
            numeric: {},
            strings: {},
        };

        const root = document.getElementById(id);
        if (root === null) {
            throw Error(`could not find HTML element #${id}`)
        }
        this._root = root;

        this._extractProperties(data.data);

        this._createPlot();
    }

    public onClick(callback: (index: number) => void) {
        this._clicked_cb = callback;
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

    private _createPlot() {
        this._plot = document.createElement("div") as unknown as PlotlyHTMLElement;
        this._plot.style.width = "100%";
        this._plot.style.height = "100%";
        this._plot.style.minHeight = "550px";
        this._root.appendChild(this._plot);

        const fullData = {
            x: this._data.numeric[this._current.x],
            y: this._data.numeric[this._current.y],
            // hovertemplate: this._color_names[0] + ": %{marker.color:.2f} <extra></extra>",
            marker: this._create_markers(),
            mode: "markers",
            type: "scattergl",
        };

        // Create a second trace to store the last clicked point, in order to
        // display it on top of the main plot with different styling
        const clicked = {
            x: [fullData.x[0]],
            y: [fullData.y[0]],
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
        }

        Plotly.newPlot(
            this._plot, [fullData as Data, clicked as Data],
            layout as Layout,
            DEFAULT_CONFIG as Config,
        );
        this._plot.on("plotly_click", (event) => this._on_plotly_click(event));
    }

    private _create_markers() {
        const rgbaColorscale: Array<[number, string]> = Inferno.map((c) => {
            return [c[0], `rgba(${c[1][0]}, ${c[1][1]}, ${c[1][2]}, 0.75)`];
        });
        const rgbColorscale: Array<[number, string]> = Inferno.map((c) => {
            return [c[0], `rgb(${c[1][0]}, ${c[1][1]}, ${c[1][2]})`];
        });

        const color = this._current.color === undefined ? 0.5 : this._data.numeric[this._current.color];
        return {
            color: color,
            colorscale: rgbaColorscale,
            line: {
                color: color,
                colorscale: rgbColorscale,
                width: 1.5,
            },
            size: 10,
            showscale: true,
            colorbar: {
                thickness: 20,
            }
        };
    }

    private _on_plotly_click(event: Plotly.PlotMouseEvent) {
        Plotly.restyle(this._plot, {
            x: [[event.points[0].x]],
            y: [[event.points[0].y]],
        }, 1);
        this._clicked_cb(event.points[0].pointNumber);
    }
}
