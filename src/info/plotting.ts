import Plotly from '../map/plotly/plotly-scatter';

import { Data } from '../map/plotly/plotly-scatter';

/**
 * wrapper around line plots of Plotly.js for 2D properties
 */
export function plotMultiDimensionalProperties(
    x: number[],
    y: number[],
    root: HTMLElement,
    plotWidth: number,
    isStatic: boolean,
    xlabel?: string,
    ylabel?: string
): void {
    const trace = {
        x: x,
        y: y,
        type: 'scatter',
        mode: 'lines',
    };

    const layout = {
        xref: 'paper',
        yref: 'paper',
        xaxis: {
            title: xlabel,
            titlefont: {
                size: 12,
            },
            showgrid: false,
            zeroline: false,
            showline: true,
            nticks: 5,
        },
        yaxis: {
            title: ylabel,
            titlefont: {
                size: 12,
            },
            showgrid: false,
            showline: true,
            zeroline: false,
            nticks: 4,
        },
        showlegend: false,
        x: 0.2,
        legend: {
            y: 0.5,
        },
        autosize: true,
        margin: {
            l: plotWidth / 6.0,
            r: 0,
            b: plotWidth / 6.0,
            t: 0,
            pad: 0,
        },
        width: plotWidth,
        height: plotWidth / 1.35, // hardcoded ratio. TODO: to be changed
        tracetoggle: false,
    };

    const config = {
        displayModeBar: false,
        responsive: true,
        staticPlot: isStatic,
    };

    const data = [trace] as Data[];

    void Plotly.newPlot(root, data, layout, config);
}
