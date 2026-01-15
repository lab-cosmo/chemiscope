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
        camera: {
            eye: { x: number; y: number; z: number };
            center: { x: number; y: number; z: number };
            up: { x: number; y: number; z: number };
            projection: { type: string };
        };
        aspectratio: { x: number; y: number; z: number };
    };
    width: number;
    height: number;
}

export interface PlotlyScatterElement extends PlotlyHTMLElement {
    _fullData: {
        type: 'scatter3d' | 'scattergl';
        uid: string;
    }[];
    _fullLayout: ScatterLayout;
}
