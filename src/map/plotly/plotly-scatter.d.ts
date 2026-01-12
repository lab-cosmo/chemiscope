export * from 'plotly.js';
import { LayoutAxis, PlotlyHTMLElement } from 'plotly.js';

interface ScatterAxis extends LayoutAxis {
    l2p(linear: number): number;
    _offset: number;
}

interface SceneAxis extends LayoutAxis {
    autorange: boolean;
}

interface ScatterLayout {
    _modeBar: {
        _uid: string;
    };

    xaxis: ScatterAxis;
    yaxis: ScatterAxis;
    scene: {
        xaxis: SceneAxis;
        yaxis: SceneAxis;
        zaxis: SceneAxis;
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
