export * from 'plotly.js';
import { LayoutAxis, PlotlyHTMLElement } from 'plotly.js';

interface ScatterAxis extends LayoutAxis {
    l2p(linear: number): number;
    _offset: number;
}

interface ScatterLayout {
    xaxis: ScatterAxis;
    yaxis: ScatterAxis;
    scene: {
        xaxis: LayoutAxis;
        yaxis: LayoutAxis;
        zaxis: LayoutAxis;
    };
}

export interface PlotlyScatterElement extends PlotlyHTMLElement {
    _fullLayout: ScatterLayout;
}
