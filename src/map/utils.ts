/**
 * @packageDocumentation
 * @module map
 */

import Plotly from './plotly/plotly-scatter';
import { PlotlyScatterElement } from './plotly/plotly-scatter';
import PNG_SVG from '../static/download-png.svg';
import SVG_SVG from '../static/download-svg.svg';
import { Warnings } from '../utils';

/**
 * Export the plot as a PNG or SVG image, hiding the "selected" trace (index 1)
 * during the export.
 *
 * @param gd Plotly plot element
 * @param format format of the image
 */
export function exportImage(gd: PlotlyScatterElement, format: 'png' | 'svg') {
    const width = Math.max(gd._fullLayout.width, 600);
    const ratio = gd._fullLayout.height / gd._fullLayout.width;

    const opts: Plotly.DownloadImgopts = {
        filename: 'chemiscope-map',
        format: format,
        width: width,
        height: format === 'png' ? width * ratio : Math.max(gd._fullLayout.height, 600),
    };

    if (format === 'png') {
        // scale is not part of `DownloadImgopts`, but accepted
        // by the function anyway
        (opts as unknown as { scale: number }).scale = 3;
    }

    // Hide the "selected" trace (index 1) for the export.
    // In 2D mode, this trace is already empty (using NaNs), but in 3D mode
    // it contains the markers for selected environments.
    Plotly.restyle(gd, { visible: false }, [1])
        .then(() => {
            return Plotly.downloadImage(gd, opts);
        })
        .then(() => {
            return Plotly.restyle(gd, { visible: true }, [1]);
        })
        .catch((e: unknown) => {
            // make sure we show the trace again even if download failed
            void Plotly.restyle(gd, { visible: true }, [1]);
            setTimeout(() => {
                throw e;
            });
        });
}

/** Extract the data associated with the first `path` element in an SVG string */
export function extractSvgPath(svg: string) {
    const doc = document.createElement('div');
    doc.innerHTML = svg;
    return doc.getElementsByTagName('path')[0].getAttribute('d');
}

/**
 * Validate min/max options provided by the user. We use `NaN` internally to mark missing values,
 * which are then transformed into undefined by this function
 */
export function getAxisRange(
    min: number,
    max: number,
    axisName: string,
    warnings: Warnings
): [number | undefined, number | undefined] {
    const minProvided = !isNaN(min);
    const maxProvided = !isNaN(max);

    // At least one range value is specified. By default, zeros are set
    if (minProvided && maxProvided) {
        if (min <= max) {
            return [min, max];
        }
        warnings.sendMessage(
            `The inserted min and max values in ${axisName} are such that min > max!` +
                `The default values will be used.`
        );
    }
    return [minProvided ? min : undefined, maxProvided ? max : undefined];
}

export function getAxisAutoRange(min: number, max: number): boolean {
    const minProvided = !isNaN(min);
    const maxProvided = !isNaN(max);

    // Both ranges are provided
    return !(minProvided && maxProvided);
}

export const DEFAULT_LAYOUT = {
    // coloraxis is used for the markers
    coloraxis: {
        cmax: 0,
        cmin: 0,
        colorbar: {
            len: 1,
            thickness: 20,
            title: {
                text: '',
                side: 'right',
                font: {
                    size: 15,
                },
            },
            y: 0,
            yanchor: 'bottom',
        },
        colorscale: [] as Plotly.ColorScale,
        showscale: true,
    },
    hovermode: 'closest',
    dragmode: 'zoom',
    legend: {
        itemclick: false,
        itemdoubleclick: false,
        tracegroupgap: 5,
        y: 1,
        yanchor: 'top',
    },
    margin: {
        b: 50,
        l: 50,
        r: 50,
        t: 50,
    },
    scene: {
        camera: {
            projection: {
                type: 'orthographic',
            },
        },
        xaxis: {
            showspikes: false,
            autorange: true,
            range: undefined as (number | undefined)[] | undefined,
            title: { text: '' },
            type: 'linear',
        },
        yaxis: {
            showspikes: false,
            autorange: true,
            range: undefined as (number | undefined)[] | undefined,
            title: { text: '' },
            type: 'linear',
        },
        zaxis: {
            showspikes: false,
            autorange: true,
            range: undefined as (number | undefined)[] | undefined,
            title: { text: '' as undefined | string },
            type: 'linear',
        },
    },
    showlegend: true,
    xaxis: {
        range: undefined as (number | undefined)[] | undefined,
        title: { text: '' },
        type: 'linear',
        zeroline: false,
    },
    yaxis: {
        range: undefined as (number | undefined)[] | undefined,
        title: { text: '' },
        type: 'linear',
        zeroline: false,
    },
    zaxis: {
        range: undefined as (number | undefined)[] | undefined,
        title: { text: '' },
        type: 'linear',
        zeroline: false,
    },
};

export const DEFAULT_CONFIG = {
    displayModeBar: true,
    displaylogo: false,
    responsive: true,
    doubleClick: 'reset',
    doubleClickDelay: 600,
    scrollZoom: true,

    modeBarButtonsToRemove: [
        'hoverClosestCartesian',
        'hoverCompareCartesian',
        'toggleSpikelines',
        'autoScale2d',
        'zoomIn2d',
        'zoomOut2d',
        'select2d',
        'lasso2d',
        'hoverClosest3d',
        'zoom3d',
        'tableRotation',
        'resetCameraLastSave3d',
        'toImage',
    ],

    modeBarButtonsToAdd: [
        [
            {
                name: 'Download PNG',
                icon: {
                    width: 400,
                    height: 447,
                    path: extractSvgPath(PNG_SVG),
                },
                click: function (gd: PlotlyScatterElement) {
                    exportImage(gd, 'png');
                },
            },
        ],
        [
            {
                name: 'Download SVG',
                icon: {
                    width: 400,
                    height: 447,
                    path: extractSvgPath(SVG_SVG),
                },
                click: function (gd: PlotlyScatterElement) {
                    exportImage(gd, 'svg');
                },
            },
        ],
    ],
};
