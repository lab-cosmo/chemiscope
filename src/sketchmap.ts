import * as Plotly from "plotly.js";
import {Config, Data, Layout, PlotlyHTMLElement} from "plotly.js";
import {Inferno} from "./colorscales";

interface SketchmapData {
    x: number[];
    y: number[];
    color_by: {
        [s: string]: number[];
    }
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
    private _name: string;
    private _data: SketchmapData;
    private _root!: HTMLElement;
    private _plot!: PlotlyHTMLElement;
    private _clicked_cb: (index: number) => void;
    private _color_names: string[];

    constructor(name: string, data: SketchmapData) {
        this._name = name;
        this._data = data;
        this._clicked_cb = (_) => { return; };
        this._color_names = [];

        for (name in this._data.color_by) {
            this._color_names.push(name);
        }
    }

    public setup(root: HTMLElement) {
        this._root = root;

        const div = document.createElement("div");
        root.appendChild(div)

        const colors = document.createElement("select") as HTMLSelectElement;
        div.appendChild(colors);
        for (let i = 0; i < this._color_names.length; i++) {
            const option = document.createElement("option") as HTMLOptionElement;
            option.value = this._color_names[i];
            option.text = this._color_names[i];
            colors.appendChild(option);
        }

        colors.onchange = () => {
            const color = this._data.color_by[colors.value];
            Plotly.restyle(this._plot, {
                'hovertemplate': colors.value + ": %{marker.color:.2f} <extra></extra>",
                'marker.showscale': true,
                'marker.color': [color],
                'marker.line': {
                    color: [color]
                },
            }, 0);
        }

        this._createPlot();
    }

    public onClick(callback: (index: number) => void) {
        this._clicked_cb = callback;
    }

    private _createPlot() {
        this._plot = document.createElement("div") as unknown as PlotlyHTMLElement;
        this._plot.style.width = "100%";
        this._plot.style.height = "100%";
        this._plot.style.minHeight = "550px";
        this._root.appendChild(this._plot);

        const fullData = {...this._data,
            hovertemplate: this._color_names[0] + ": %{marker.color:.2f} <extra></extra>",
            marker: this._create_markers(),
            mode: "markers",
            type: "scattergl",
        };

        // Create a second trace to store the last clicked point, in order to
        // display it on top of the main plot with different styling
        const clicked = {
            x: [this._data.x[0]],
            y: [this._data.y[0]],
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

        const color = this._data.color_by[this._color_names[0]];
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
