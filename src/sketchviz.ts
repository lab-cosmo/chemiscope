import Plotly from "plotly.js";
import {Config, Data, Layout, PlotlyHTMLElement} from "plotly.js";
import {Inferno} from "./colorscales";

interface SketchvizData {
    x: number[];
    y: number[];
    color: number[];
}

export class Sketchviz {
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

    private static create_markers(colors: number[]) {
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

    private name: string;
    private data: SketchvizData;
    // private structure: {};

    private root: HTMLElement;
    private plot: PlotlyHTMLElement;

    constructor(name: string, data: SketchvizData, _: {}) {
        this.name = name;
        this.data = data;
        // this.structure = structure;
    }

    public setup(root: HTMLElement) {
        this.root = root;

        const title = document.createElement("h3");
        title.innerHTML = this.name;
        this.root.appendChild(title);

        this.createPlot();
    }

    private createPlot() {
        this.plot = document.createElement("div") as unknown as PlotlyHTMLElement;
        this.root.appendChild(this.plot);

        this.plot.setAttribute("style", "width: 600px; height: 600px;");

        const fullData = {...this.data,
            hoverinfo: "none",
            marker: Sketchviz.create_markers(this.data.color),
            mode: "markers",
            type: "scattergl",
        };

        // Create a second trace to store the last clicked point, in order to
        // display it on top of the main plot with different styling
        const clicked = {
            x: [NaN],
            y: [NaN],
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
            this.plot, [fullData as Data, clicked as Data],
            Sketchviz.LAYOUT as Layout,
            Sketchviz.CONFIG as Config,
        );
        this.plot.on("plotly_click", (event) => this.on_plotly_click(event));
    }

    private on_plotly_click(data: any) {
        Plotly.restyle(this.plot, {
            x: [[data.points[0].x]],
            y: [[data.points[0].y]],
        }, 1);
    }
}
