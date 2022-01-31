export * from 'plotly.js';
import { LayoutAxis, PlotlyHTMLElement } from 'plotly.js';

interface ScatterAxis extends LayoutAxis {
    l2p(linear: number): number;
    _offset: number;
}

interface ScatterLayout {
    _modeBar: {
        _uid: string;
    };

    xaxis: ScatterAxis;
    yaxis: ScatterAxis;
    scene: {
        xaxis: LayoutAxis;
        yaxis: LayoutAxis;
        zaxis: LayoutAxis;
    };
    width: number;
    height: number;
}

export interface PlotlyScatterElement extends PlotlyHTMLElement {
    _fullData: {
        type: 'scatter3d' | 'scattergl';
    }[];
    _fullLayout: ScatterLayout;
}
