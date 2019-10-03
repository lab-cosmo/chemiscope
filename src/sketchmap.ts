import * as Plotly from "plotly.js";
import {Config, Data, Layout, PlotlyHTMLElement} from "plotly.js";
import {Inferno} from "./colorscales";

interface SketchmapData {
    x: number[];
    y: number[];
    color: number[];
}

export class Sketchmap {
    private static LAYOUT = {
        hovermode: "closest",
        showlegend: false,
    };

    private static CONFIG = {
        displayModeBar: true,
        displaylogo: false,
        modeBarButtonsToRemove: [
            "hoverClosestCartesian",
            "hoverCompareCartesian",
            "toggleSpikelines",
        ],
    };

    private static _create_markers(colors: number[]) {
        const rgbaColorscale: Array<[number, string]> = Inferno.map((c) => {
            return [c[0], `rgba(${c[1][0]}, ${c[1][1]}, ${c[1][2]}, 0.75)`];
        });
        const rgbColorscale: Array<[number, string]> = Inferno.map((c) => {
            return [c[0], `rgb(${c[1][0]}, ${c[1][1]}, ${c[1][2]})`];
        });

        return {
            color: colors,
            colorscale: rgbaColorscale,
            line: {
                color: colors,
                colorscale: rgbColorscale,
                width: 1.5,
            },
            size: 10,
        };
    }

    private _name: string;
    private _data: SketchmapData;
    private _root: HTMLElement;
    private _plot: PlotlyHTMLElement;
    private _clicked_cb: (index: number) => void;

    constructor(name: string, data: SketchmapData) {
        this._name = name;
        this._data = data;
        this._clicked_cb = (_) => { return; };
    }

    public setup(root: HTMLElement) {
        this._root = root;

        const title = document.createElement("h3");
        title.innerHTML = this._name;
        this._root.appendChild(title);

        this._createPlot();
    }

    public onClick(callback: (index: number) => void) {
        this._clicked_cb = callback;
    }

    private _createPlot() {
        this._plot = document.createElement("div") as unknown as PlotlyHTMLElement;
        this._root.appendChild(this._plot);

        this._plot.setAttribute("style", "width: 500px; height: 500px;");

        const fullData = {...this._data,
            hoverinfo: "none",
            marker: Sketchmap._create_markers(this._data.color),
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
                color: "blue",
                line: {
                    color: "black",
                    width: 1.5,
                },
                size: 18,
            },
            hoverinfo: "none",
        };

        Plotly.newPlot(
            this._plot, [fullData as Data, clicked as Data],
            Sketchmap.LAYOUT as Layout,
            Sketchmap.CONFIG as Config,
        );
        this._plot.on("plotly_click", (event) => this._on_plotly_click(event));
    }

    private _on_plotly_click(data: Plotly.PlotMouseEvent) {
        Plotly.restyle(this._plot, {
            x: [[data.points[0].x]],
            y: [[data.points[0].y]],
        }, 1);
        this._clicked_cb(data.points[0].pointNumber);
    }
}
