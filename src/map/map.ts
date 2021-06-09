/**
 * @packageDocumentation
 * @module map
 */

import assert from 'assert';

import Plotly from './plotly/plotly-scatter';
import { Config, Data, Layout, PlotlyScatterElement } from './plotly/plotly-scatter';

import { Property } from '../dataset';

import { EnvironmentIndexer, Indexes } from '../indexer';
import { OptionModificationOrigin, SavedSettings } from '../options';
import { GUID, PositioningCallback, arrayMaxMin } from '../utils';
import { enumerate, getByID, getFirstKey } from '../utils';

import { MapData, NumericProperty } from './data';
import { MarkerData } from './marker';
import { AxisOptions, MapOptions, get3DSymbol } from './options';

// var extract = require('extract-svg-path');
// var parse = require('parse-svg-path');
// import PNG_SVG from '../static/png-icon.svg';
// var path_PNG = extract(PNG_SVG);
// var svg = parse(path_PNG);
// import SVG_SVG from '../static/svg-icon.svg';

const DEFAULT_LAYOUT = {
    // coloraxis is used for the markers
    coloraxis: {
        cmax: 0,
        cmin: 0,
        colorbar: {
            len: 1,
            thickness: 20,
            title: {
                text: '' as undefined | string,
            },
            y: 0,
            yanchor: 'bottom',
        },
        colorscale: [] as Plotly.ColorScale,
        showscale: true,
    },
    hovermode: 'closest',
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
            title: '',
        },
        yaxis: {
            showspikes: false,
            title: '',
        },
        zaxis: {
            showspikes: false,
            title: '' as undefined | string,
        },
    },
    showlegend: true,
    xaxis: {
        range: undefined,
        title: '',
        type: 'linear',
        zeroline: false,
    },
    yaxis: {
        range: undefined,
        title: '',
        type: 'linear',
        zeroline: false,
    },
};

const DEFAULT_CONFIG = {
    displayModeBar: true,
    displaylogo: false,
    responsive: true,
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
        'tableRotation',
        'resetCameraLastSave3d',
        'toImage',
    ],

    modeBarButtonsToAdd: [
        {
            name: 'Download as PNG',
            icon: {
                width: 550,
                height: 550,
                path: 'M488.426,197.019H475.2v-63.816c0-0.398-0.063-0.799-0.116-1.202c-0.021-2.534-0.827-5.023-2.562-6.995L366.325,3.694 c-0.032-0.031-0.063-0.042-0.085-0.076c-0.633-0.707-1.371-1.295-2.151-1.804c-0.231-0.155-0.464-0.285-0.706-0.419 c-0.676-0.369-1.393-0.675-2.131-0.896c-0.2-0.056-0.38-0.138-0.58-0.19C359.87,0.119,359.037,0,358.193,0H97.2 c-11.918,0-21.6,9.693-21.6,21.601v175.413H62.377c-17.049,0-30.873,13.818-30.873,30.873v160.545 c0,17.043,13.824,30.87,30.873,30.87h13.224V529.2c0,11.907,9.682,21.601,21.6,21.601h356.4c11.907,0,21.6-9.693,21.6-21.601 V419.302h13.226c17.044,0,30.871-13.827,30.871-30.87v-160.54C519.297,210.838,505.47,197.019,488.426,197.019z M97.2,21.605 h250.193v110.513c0,5.967,4.841,10.8,10.8,10.8h95.407v54.108H97.2V21.605z M338.871,225.672L284.545,386.96h-42.591 l-51.69-161.288h39.967l19.617,68.196c5.508,19.143,10.531,37.567,14.36,57.67h0.717c4.061-19.385,9.089-38.527,14.592-56.953 l20.585-68.918h38.77V225.672z M68.458,379.54l7.415-30.153c9.811,5.021,24.888,10.051,40.439,10.051 c16.751,0,25.607-6.935,25.607-17.465c0-10.052-7.662-15.795-27.05-22.734c-26.8-9.328-44.263-24.168-44.263-47.611 c0-27.524,22.971-48.579,61.014-48.579c18.188,0,31.591,3.823,41.159,8.131l-8.126,29.437c-6.465-3.116-17.945-7.657-33.745-7.657 c-15.791,0-23.454,7.183-23.454,15.552c0,10.296,9.089,14.842,29.917,22.731c28.468,10.536,41.871,25.365,41.871,48.094 c0,27.042-20.812,50.013-65.09,50.013C95.731,389.349,77.538,384.571,68.458,379.54z M453.601,523.353H97.2V419.302h356.4V523.353z M488.911,379.54c-11.243,3.823-32.537,9.103-53.831,9.103c-29.437,0-50.73-7.426-65.57-21.779 c-14.839-13.875-22.971-34.942-22.738-58.625c0.253-53.604,39.255-84.235,92.137-84.235c20.81,0,36.852,4.073,44.74,7.896 l-7.657,29.202c-8.859-3.829-19.849-6.95-37.567-6.95c-30.396,0-53.357,17.233-53.357,52.173c0,33.265,20.81,52.882,50.73,52.882 c8.375,0,15.072-0.96,17.94-2.395v-33.745h-24.875v-28.471h60.049V379.54L488.911,379.54z',
            },
            click: function (gd: PlotlyScatterElement) {
                Plotly.downloadImage(gd, {
                    filename: 'newplot',
                    format: 'png',
                    width: gd._fullLayout.width,
                    height: gd._fullLayout.height,
                });
            },
        },
        {
            name: 'Download as SVG',
            icon: {
                width: 550,
                height: 550,
                path: 'M488.426,197.019H475.2v-63.816c0-0.398-0.063-0.799-0.116-1.202c-0.021-2.534-0.827-5.023-2.562-6.995L366.325,3.694 c-0.032-0.031-0.063-0.042-0.085-0.076c-0.633-0.707-1.371-1.295-2.151-1.804c-0.231-0.155-0.464-0.285-0.706-0.419 c-0.676-0.369-1.393-0.675-2.131-0.896c-0.2-0.056-0.38-0.138-0.58-0.19C359.87,0.119,359.037,0,358.193,0H97.2 c-11.918,0-21.6,9.693-21.6,21.601v175.413H62.377c-17.049,0-30.873,13.818-30.873,30.873v160.545 c0,17.043,13.824,30.87,30.873,30.87h13.224V529.2c0,11.907,9.682,21.601,21.6,21.601h356.4c11.907,0,21.6-9.693,21.6-21.601 V419.302h13.226c17.044,0,30.871-13.827,30.871-30.87v-160.54C519.297,210.838,505.47,197.019,488.426,197.019z M97.2,21.605 h250.193v110.513c0,5.967,4.841,10.8,10.8,10.8h95.407v54.108H97.2V21.605z M338.871,225.672L284.545,386.96h-42.591 l-51.69-161.288h39.967l19.617,68.196c5.508,19.143,10.531,37.567,14.36,57.67h0.717c4.061-19.385,9.089-38.527,14.592-56.953 l20.585-68.918h38.77V225.672z M68.458,379.54l7.415-30.153c9.811,5.021,24.888,10.051,40.439,10.051 c16.751,0,25.607-6.935,25.607-17.465c0-10.052-7.662-15.795-27.05-22.734c-26.8-9.328-44.263-24.168-44.263-47.611 c0-27.524,22.971-48.579,61.014-48.579c18.188,0,31.591,3.823,41.159,8.131l-8.126,29.437c-6.465-3.116-17.945-7.657-33.745-7.657 c-15.791,0-23.454,7.183-23.454,15.552c0,10.296,9.089,14.842,29.917,22.731c28.468,10.536,41.871,25.365,41.871,48.094 c0,27.042-20.812,50.013-65.09,50.013C95.731,389.349,77.538,384.571,68.458,379.54z M453.601,523.353H97.2V419.302h356.4V523.353z M488.911,379.54c-11.243,3.823-32.537,9.103-53.831,9.103c-29.437,0-50.73-7.426-65.57-21.779 c-14.839-13.875-22.971-34.942-22.738-58.625c0.253-53.604,39.255-84.235,92.137-84.235c20.81,0,36.852,4.073,44.74,7.896 l-7.657,29.202c-8.859-3.829-19.849-6.95-37.567-6.95c-30.396,0-53.357,17.233-53.357,52.173c0,33.265,20.81,52.882,50.73,52.882 c8.375,0,15.072-0.96,17.94-2.395v-33.745h-24.875v-28.471h60.049V379.54L488.911,379.54z',
            },
            click: function (gd: PlotlyScatterElement) {
                Plotly.downloadImage(gd, {
                    filename: 'newplot',
                    format: 'svg',
                    width: gd._fullLayout.width,
                    height: gd._fullLayout.height,
                });
            },
        },
    ],
};

/**
 * The [[PropertiesMap]] class displays a 2D or 3D map (scatter plot) of
 * properties in the dataset, using [plotly.js](https://plot.ly/javascript/)
 * for rendering.
 *
 * Properties can be used as x, y, or z values, as well as points color and
 * size. Additionally, string properties can be used as symbols for the scatter
 * plot markers.
 */
export class PropertiesMap {
    /** Callback fired when the plot is clicked and the position of the active marker changes */
    public onselect: (indexes: Indexes) => void;
    /**
     * Callback fired when the active marker is changed by clicking on the map
     *
     * @param guid GUID of the new active marker
     * @param indexes index of the environment the new active marker is showing
     */
    public activeChanged: (guid: GUID, indexes: Indexes) => void;

    /**
     * Callback to get the initial positioning of the settings modal.
     *
     * The callback gets the current placement of the settings as a
     * [DOMRect](https://developer.mozilla.org/en-US/docs/Web/API/DOMRect), and
     * should return top and left positions in pixels, used with `position:
     * fixed`. The callback is called once, the first time the settings are
     * opened.
     */
    public positionSettingsModal: PositioningCallback;

    /// HTML root holding the full plot
    private _root: HTMLElement;
    /// Plotly plot
    private _plot!: PlotlyScatterElement;
    /// All known properties
    private _data: MapData;

    /// GUID of the currently selected point
    private _active?: GUID;
    /// Map associating currently selected markers GUID to additional data
    private _selected: Map<GUID, MarkerData>;

    /// environment indexer
    private _indexer: EnvironmentIndexer;
    /// Settings of the map
    private _options: MapOptions;
    /// Button used to reset the range of color axis
    private _colorReset: HTMLButtonElement;

    /**
     * Create a new [[PropertiesMap]] inside the DOM element with the given HTML
     * `id`
     *
     * @param id         HTML id of the DOM element where the map should live
     * @param indexer    [[EnvironmentIndexer]] used to translate indexes from
     *                   environments index to structure/atom indexes
     * @param properties properties to be displayed
     */
    constructor(
        config: { id: string; settings: SavedSettings },
        indexer: EnvironmentIndexer,
        properties: { [name: string]: Property }
    ) {
        this._indexer = indexer;
        this.onselect = () => {};
        this.activeChanged = () => {};
        this._selected = new Map<GUID, MarkerData>();

        this._root = getByID(config.id);
        if (this._root.style.position === '') {
            this._root.style.position = 'relative';
        }

        this._plot = document.createElement('div') as unknown as PlotlyScatterElement;
        this._plot.style.width = '100%';
        this._plot.style.height = '100%';
        this._root.appendChild(this._plot);

        this._data = new MapData(properties);

        this._options = new MapOptions(
            this._root,
            this._data[this._indexer.mode],
            (rect) => this.positionSettingsModal(rect),
            config.settings
        );
        this._colorReset = getByID<HTMLButtonElement>('chsp-color-reset');

        this._connectSettings();

        // By default, position the modal for settings on top of the plot,
        // centered horizontally
        this.positionSettingsModal = (rect: DOMRect) => {
            const rootRect = this._root.getBoundingClientRect();
            return {
                left: rootRect.left + rootRect.width / 2 - rect.width / 2,
                top: rootRect.top + 20,
            };
        };

        this._createPlot();
    }

    /**
     * Remove all HTML added by this [[PropertiesMap]] in the current document
     */
    public remove(): void {
        this._root.innerHTML = '';
        this._options.remove();
    }

    /**
     * Change the environment indicated by the currently active marker to
     * the one with the given `indexes`
     */
    public select(indexes: Indexes): void {
        if (this._active === undefined) {
            throw Error('tries to update selected environment, but there is no active marker');
        }

        const data = this._selected.get(this._active);
        assert(data !== undefined);
        // this prevents infinite recursion
        // will return false if the data already corresponds to this index
        if (data.select(indexes)) {
            this._updateMarkers();
        }
    }

    /**
     * Set the marker with given GUID as the active marker.
     *
     * @param guid the GUID of the new active viewer
     */
    public setActive(guid: GUID): void {
        if (this._active !== undefined) {
            const oldData = this._selected.get(this._active);
            assert(oldData !== undefined);
            oldData.deactivate();
        }

        this._active = guid;
        const data = this._selected.get(this._active);
        assert(data !== undefined);
        data.activate();

        if (this._is3D()) {
            this._restyle({ 'marker.size': this._sizes(1) } as Data, 1);
        }
    }

    /**
     * Add a new marker to the map. The new marker is set as the active one
     *
     * @param guid GUID of the new marker
     * @param indexes indexes of the environment that the new marker should show
     */
    public addMarker(guid: GUID, color: string, indexes: Indexes): void {
        assert(!this._selected.has(guid));

        const data = new MarkerData(guid, color, indexes.environment);
        this._root.appendChild(data.marker);
        data.marker.onclick = () => {
            this.setActive(guid);
            this.activeChanged(guid, this._indexer.from_environment(data.current));
        };
        this._selected.set(guid, data);
        this._updateMarkers([data]);
        this.setActive(guid);
    }

    /**
     * Removes a marker from the map.
     *
     * @param guid GUID of the marker to remove
     */
    public removeMarker(guid: GUID): void {
        // If we remove the active marker, let's change which one is active
        if (this._active === guid) {
            if (this._selected.size === 1) {
                // we are removing the last one
                this._active = undefined;
            } else {
                this.setActive(getFirstKey(this._selected, guid));
            }
        }

        // remove HTML marker
        const data = this._selected.get(guid);
        assert(data !== undefined);
        data.remove();

        this._selected.delete(guid);
        this._updateMarkers();
    }

    /**
     * Apply saved settings to the map.
     */
    public applySettings(settings: SavedSettings): void {
        this._options.applySettings(settings);
    }

    /**
     * Save the values of the current settings in a way that an be used with
     * [[applySettings]] or saved to JSON.
     */
    public saveSettings(): SavedSettings {
        return this._options.saveSettings();
    }

    /** Forward to Ploty.restyle */
    private _restyle(data: Partial<Data>, traces?: number | number[]) {
        Plotly.restyle(this._plot, data, traces).catch((e) =>
            setTimeout(() => {
                throw e;
            })
        );
    }

    /** Forward to Ploty.relayout */
    private _relayout(layout: Partial<Layout>) {
        Plotly.relayout(this._plot, layout).catch((e) =>
            setTimeout(() => {
                throw e;
            })
        );
    }

    /** Add all the required callback to the settings */
    private _connectSettings() {
        // ======= x axis settings
        this._options.x.property.onchange = () => {
            const values = this._coordinates(this._options.x) as number[][];
            this._restyle({ x: values }, [0, 1]);
            this._relayout({
                'scene.xaxis.title': this._title(this._options.x.property.value),
                'xaxis.title': this._title(this._options.x.property.value),
            } as unknown as Layout);
        };

        this._options.x.scale.onchange = () => {
            if (this._is3D()) {
                this._relayout({
                    'scene.xaxis.type': this._options.x.scale.value,
                } as unknown as Layout);
            } else {
                this._relayout({ 'xaxis.type': this._options.x.scale.value as Plotly.AxisType });
            }
        };

        // function creating a function to be used as onchange callback
        // for <axis>.min and <axis>.max
        const rangeChange = (name: string, axis: AxisOptions) => {
            return (_: number, origin: OptionModificationOrigin) => {
                if (origin === 'JS') {
                    // prevent recursion: this function calls relayout, which then
                    // calls _afterplot, which reset the min/max values.
                    return;
                }
                const min = axis.min.value;
                const max = axis.max.value;
                if (this._is3D()) {
                    this._relayout({
                        [`scene.${name}.range`]: [min, max],
                    } as unknown as Layout);
                } else {
                    this._relayout({ [`${name}.range`]: [min, max] });
                }
            };
        };

        this._options.x.min.onchange = rangeChange('xaxis', this._options.x);
        this._options.x.max.onchange = rangeChange('xaxis', this._options.x);

        // ======= y axis settings
        this._options.y.property.onchange = () => {
            const values = this._coordinates(this._options.y) as number[][];
            this._restyle({ y: values }, [0, 1]);
            this._relayout({
                'scene.yaxis.title': this._title(this._options.y.property.value),
                'yaxis.title': this._title(this._options.y.property.value),
            } as unknown as Layout);
        };

        this._options.y.scale.onchange = () => {
            if (this._is3D()) {
                this._relayout({
                    'scene.yaxis.type': this._options.y.scale.value,
                } as unknown as Layout);
            } else {
                this._relayout({ 'yaxis.type': this._options.y.scale.value as Plotly.AxisType });
            }
        };

        this._options.y.min.onchange = rangeChange('yaxis', this._options.y);
        this._options.y.max.onchange = rangeChange('yaxis', this._options.y);

        // ======= z axis settings
        this._options.z.property.onchange = () => {
            // eslint-disable-next-line @typescript-eslint/no-unsafe-member-access, @typescript-eslint/no-explicit-any
            const was3D = (this._plot as any)._fullData[0].type === 'scatter3d';
            if (this._options.z.property.value === '') {
                if (was3D) {
                    this._switch2D();
                }
            } else {
                if (!was3D) {
                    this._switch3D();
                }
            }

            const values = this._coordinates(this._options.z);
            this._restyle({ z: values } as Data, [0, 1]);
            this._relayout({
                'scene.zaxis.title': this._title(this._options.z.property.value),
            } as unknown as Layout);
        };

        this._options.z.scale.onchange = () => {
            if (this._options.z.property.value !== '') {
                this._relayout({
                    'scene.zaxis.type': this._options.z.scale.value,
                } as unknown as Layout);
            }
        };

        this._options.z.min.onchange = rangeChange('zaxis', this._options.z);
        this._options.z.max.onchange = rangeChange('zaxis', this._options.z);

        // ======= color axis settings
        // setup initial state of the color settings
        if (this._options.color.property.value) {
            this._options.color.enable();
            this._colorReset.disabled = false;

            const values = this._colors(0)[0] as number[];
            const { min, max } = arrayMaxMin(values);

            this._options.color.min.value = min;
            this._options.color.max.value = max;
        } else {
            this._options.color.min.disable();
            this._colorReset.disabled = true;

            this._options.color.min.value = 0;
            this._options.color.max.value = 0;
        }

        this._options.color.property.onchange = () => {
            if (this._options.color.property.value !== '') {
                this._options.color.enable();
                this._colorReset.disabled = false;

                const values = this._colors(0)[0] as number[];
                const { min, max } = arrayMaxMin(values);

                this._options.color.min.value = min;
                this._options.color.max.value = max;

                this._relayout({
                    'coloraxis.colorbar.title.text': this._title(
                        this._options.color.property.value
                    ),
                    'coloraxis.showscale': true,
                } as unknown as Layout);
            } else {
                this._options.color.disable();
                this._colorReset.disabled = true;

                this._options.color.min.value = 0;
                this._options.color.max.value = 0;

                this._relayout({
                    'coloraxis.colorbar.title.text': undefined,
                    'coloraxis.showscale': false,
                } as unknown as Layout);
            }

            this._restyle(
                {
                    hovertemplate: this._options.hovertemplate(),
                    'marker.color': this._colors(0),
                } as Data,
                0
            );
        };

        const colorRangeChange = () => {
            const min = this._options.color.min.value;
            const max = this._options.color.max.value;
            this._relayout({
                'coloraxis.cmax': max,
                'coloraxis.cmin': min,
                // looks like changing only 'coloraxis.cmax'/'coloraxis.cmin' do
                // not update the color of the points (although it does change
                // the colorbar). Asking for an update of 'coloraxis.colorscale'
                // seems to do the trick. This is possiblely a Ploty bug, we
                // would need to investiguate a bit more.
                'coloraxis.colorscale': this._options.colorScale(),
            } as unknown as Layout);
        };
        this._options.color.min.onchange = colorRangeChange;
        this._options.color.max.onchange = colorRangeChange;

        this._colorReset.onclick = () => {
            const values = this._colors(0)[0] as number[];
            const { min, max } = arrayMaxMin(values);
            this._options.color.min.value = min;
            this._options.color.max.value = max;
            this._relayout({
                'coloraxis.cmax': max,
                'coloraxis.cmin': min,
                // same as above regarding update of the points color
                'coloraxis.colorscale': this._options.colorScale(),
            } as unknown as Layout);
        };

        // ======= color palette
        this._options.palette.onchange = () => {
            this._relayout({
                'coloraxis.colorscale': this._options.colorScale(),
            } as unknown as Layout);
        };

        // ======= markers symbols
        this._options.symbol.onchange = () => {
            this._restyle({ 'marker.symbol': this._symbols() }, [0, 1]);

            this._restyle({
                name: this._legendNames(),
                showlegend: this._showlegend(),
            } as unknown as Data);

            this._relayout({
                'coloraxis.colorbar.len': this._colorbarLen(),
            } as unknown as Layout);
        };

        // ======= markers size
        this._options.size.mode.onchange = () => {
            if (this._options.size.mode.value !== 'constant') {
                this._options.size.property.enable();
                this._options.size.reverse.enable();
            } else {
                this._options.size.property.disable();
                this._options.size.reverse.disable();
            }
            this._restyle({ 'marker.size': this._sizes(0) } as Data, 0);
        };

        this._options.size.factor.onchange = () => {
            this._restyle({ 'marker.size': this._sizes(0) } as Data, 0);
        };

        this._options.size.property.onchange = () => {
            this._restyle({ 'marker.size': this._sizes(0) } as Data, 0);
        };

        this._options.size.reverse.onchange = () => {
            this._restyle({ 'marker.size': this._sizes(0) } as Data, 0);
        };
    }

    /** Actually create the Plotly plot */
    private _createPlot() {
        this._plot.innerHTML = '';

        const type = this._is3D() ? 'scatter3d' : 'scattergl';
        // The main trace, containing default data
        const main = {
            name: '',
            type: type,

            x: this._coordinates(this._options.x, 0)[0],
            y: this._coordinates(this._options.y, 0)[0],
            z: this._coordinates(this._options.z, 0)[0],

            hovertemplate: this._options.hovertemplate(),
            marker: {
                color: this._colors(0)[0],
                coloraxis: 'coloraxis',
                line: {
                    color: this._lineColors(0)[0],
                    width: 1,
                },
                // prevent plolty from messing with opacity when doing bubble
                // style charts (different sizes for each point)
                opacity: 1,
                size: this._sizes(0)[0],
                sizemode: 'area',
                symbol: this._symbols(0)[0],
            },
            mode: 'markers',
            showlegend: false,
        };

        // Create a second trace to store the last clicked point, in order to
        // display it on top of the main plot with different styling. This is
        // only used in 3D mode, since it is way slower than moving
        // this._selectedMarker around.
        const selected = {
            name: 'selected',
            type: type,

            x: [],
            y: [],
            z: [],

            hoverinfo: 'none',
            marker: {
                color: [],
                line: {
                    color: [],
                    width: 2,
                },
                opacity: 1,
                size: [],
                sizemode: 'area',
            },
            mode: 'markers',
            showlegend: false,
        };

        const traces = [main as Data, selected as Data];

        const legendNames = this._legendNames().slice(2);
        const showlegend = this._showlegend().slice(2);
        assert(legendNames.length === showlegend.length);
        const currentLength = legendNames.length;

        if (this._data.maxSymbols > 0) {
            // resize & fill arrays
            legendNames.length = this._data.maxSymbols;
            legendNames.fill('', currentLength);
            showlegend.length = this._data.maxSymbols;
            showlegend.fill(false, currentLength);
        }

        // add empty traces to be able to display the symbols legend
        // one trace for each possible symbol
        for (let s = 0; s < this._data.maxSymbols; s++) {
            const data = {
                name: legendNames[s],
                type: type,

                // We need to add a dummy point to force plotly to display the
                // associated legend; but we don't want to see the point in the
                // map. Setting the coordinates to NaN acheive this.
                x: [NaN],
                y: [NaN],
                z: [NaN],

                marker: {
                    color: 'black',
                    size: 10,
                    symbol: this._is3D() ? get3DSymbol(s) : s,
                },
                mode: 'markers',
                showlegend: showlegend[s],
            };
            traces.push(data as Data);
        }

        // make a copy of the default layout
        const layout = JSON.parse(JSON.stringify(DEFAULT_LAYOUT)) as typeof DEFAULT_LAYOUT;
        // and set values speific to the displayed dataset
        layout.xaxis.title = this._title(this._options.x.property.value);
        layout.yaxis.title = this._title(this._options.y.property.value);
        layout.xaxis.type = this._options.x.scale.value;
        layout.yaxis.type = this._options.y.scale.value;
        layout.scene.xaxis.title = this._title(this._options.x.property.value);
        layout.scene.yaxis.title = this._title(this._options.y.property.value);
        layout.scene.zaxis.title = this._title(this._options.z.property.value);
        layout.coloraxis.colorscale = this._options.colorScale();
        layout.coloraxis.cmin = this._options.color.min.value;
        layout.coloraxis.cmax = this._options.color.max.value;
        layout.coloraxis.colorbar.title.text = this._title(this._options.color.property.value);
        layout.coloraxis.colorbar.len = this._colorbarLen();

        // Create an empty plot and fill it below
        Plotly.newPlot(
            this._plot,
            traces,
            layout as Partial<Layout>,
            DEFAULT_CONFIG as unknown as Config
        ).catch((e) =>
            setTimeout(() => {
                throw e;
            })
        );

        this._plot.on('plotly_click', (event: Plotly.PlotMouseEvent) => {
            // don't update selected env on double click, since it is bound to
            // 'reset zoom level' in 2D mode.
            if (event.event && event.event.detail === 2) {
                return;
            }

            let environment = event.points[0].pointNumber;
            if (this._is3D() && event.points[0].data.name === 'selected') {
                // if someone has clicked on a selection marker, set to active
                // this is only used in 3D mode, since in 2D the HTML marker
                // directly deal with the click event
                for (const [i, [guid, data]] of enumerate(this._selected.entries())) {
                    if (event.points[0].pointNumber === i) {
                        environment = data.current;
                        if (this._active !== guid) {
                            this.setActive(guid);
                            this.activeChanged(guid, this._indexer.from_environment(data.current));
                        }
                        break;
                    }
                }
            }

            const indexes = this._indexer.from_environment(environment);

            this.select(indexes);
            this.onselect(indexes);
        });

        this._plot.on('plotly_afterplot', () => this._afterplot());
        this._updateMarkers();
    }

    /** Get the property with the given name */
    private _property(name: string): NumericProperty {
        const result = this._data[this._indexer.mode][name];
        if (result === undefined) {
            throw Error(`unknown property '${name}' requested in map`);
        }
        return result;
    }

    /**
     * Get the values associated with the given `axis`, to use with the given
     * plotly `trace`, or all of them if `trace === undefined`
     *
     * @param  axis   Options of the axis we need coordinates for
     * @param  trace  plotly trace for which we require coordinate
     * @return        data usable with Plotly.restyle
     */
    private _coordinates(axis: AxisOptions, trace?: number): Array<undefined | number[]> {
        // this happen for the z axis in 2D mode
        if (axis.property.value === '') {
            return this._selectTrace(undefined, undefined, trace);
        }

        const values = this._property(axis.property.value).values;
        // in 2d mode, set all selected markers coordinates to NaN since we are
        // using HTML markers instead.
        const selected = [];
        for (const marker of this._selected.values()) {
            if (this._is3D()) {
                selected.push(values[marker.current]);
            } else {
                selected.push(NaN);
            }
        }
        return this._selectTrace<number[]>(values, selected, trace);
    }

    private _title(name: string): string {
        if (name !== '') {
            const units = this._property(name).units;
            if (units !== undefined) {
                return name + `/${units}`;
            }
        }
        return name;
    }

    /**
     * Get the color values to use with the given plotly `trace`, or all of
     * them if `trace === undefined`
     */
    private _colors(trace?: number): Array<string | string[] | number | number[]> {
        let values;
        if (this._options.hasColors()) {
            values = this._property(this._options.color.property.value).values;
        } else {
            values = 0.5;
        }

        const selected = [];
        for (const data of this._selected.values()) {
            selected.push(data.color);
        }

        return this._selectTrace<string | string[] | number | number[]>(values, selected, trace);
    }

    /**
     * Get the **line** color values to use with the given plotly `trace`, or
     * all of them if `trace === undefined`
     */
    private _lineColors(trace?: number): string[] {
        if (this._is3D()) {
            return this._selectTrace<string>('black', 'black', trace);
        } else {
            return this._selectTrace<string>('rgba(1, 1, 1, 0.3)', 'black', trace);
        }
    }

    /**
     * Get the values to use as marker size with the given plotly `trace`, or
     * all of them if `trace === undefined`.
     */
    private _sizes(trace?: number): Array<number | number[]> {
        const sizes = this._property(this._options.size.property.value).values;
        const values = this._options.calculateSizes(sizes);
        const selected = [];
        if (this._is3D()) {
            for (const guid of this._selected.keys()) {
                if (guid === this._active) {
                    selected.push(600);
                } else {
                    selected.push(300);
                }
            }
        }
        return this._selectTrace<number | number[]>(values, selected, trace);
    }

    /**
     * Get the values to use as marker symbol with the given plotly `trace`, or
     * all of them if `trace === undefined`.
     */
    private _symbols(trace?: number): Array<string | string[] | number[]> {
        if (this._options.symbol.value === '') {
            // default to 0 (i.e. circles)
            return this._selectTrace<string | string[]>('circle', 'circle', trace);
        }

        const property = this._property(this._options.symbol.value);
        const symbols = this._options.getSymbols(property);
        const selected = [];
        for (const data of this._selected.values()) {
            selected.push(symbols[data.current]);
        }
        return this._selectTrace<typeof symbols>(symbols, selected as typeof symbols, trace);
    }

    /** Should we show the legend for the various symbols used? */
    private _showlegend(): boolean[] {
        const result = [false, false];

        if (this._options.symbol.value !== '') {
            for (let i = 0; i < this._symbolsCount(); i++) {
                result.push(true);
            }
            return result;
        } else {
            for (let i = 0; i < this._data.maxSymbols; i++) {
                result.push(false);
            }
            return result;
        }
    }

    /** Get the list of symbols names to use for the legend */
    private _legendNames(): string[] {
        const result = ['', ''];

        if (this._options.symbol.value !== '') {
            const property = this._property(this._options.symbol.value);
            assert(property.string !== undefined);
            for (const name of property.string.strings()) {
                result.push(name);
            }
            return result;
        } else {
            for (let i = 0; i < this._data.maxSymbols; i++) {
                result.push('');
            }
            return result;
        }
    }

    /**
     * Select either main, selected or both depending on `trace`, and return
     * them in a mode usable with `Plotly.restyle`/[[PropertiesMap._restyle]]
     */
    private _selectTrace<T>(main: T, selected: T, trace?: number): T[] {
        if (trace === 0) {
            return [main];
        } else if (trace === 1) {
            return [selected];
        } else if (trace === undefined) {
            return [main, selected];
        } else {
            throw Error('internal error: invalid trace number');
        }
    }

    /** Get the length of the colorbar to accomodate for the legend */
    private _colorbarLen(): number {
        /// Heigh of a legend item in plot unit
        const LEGEND_ITEM_HEIGH = 0.045;
        return 1 - LEGEND_ITEM_HEIGH * this._symbolsCount();
    }

    /** Is the the current plot a 3D plot? */
    private _is3D(): boolean {
        return this._options.is3D();
    }

    /** How many symbols are on this plot?*/
    private _symbolsCount(): number {
        if (this._options.symbol.value !== '') {
            const property = this._property(this._options.symbol.value);
            assert(property.string !== undefined);
            return property.string.strings().length;
        } else {
            return 0;
        }
    }

    /** Switch current plot from 2D to 3D */
    private _switch3D(): void {
        assert(this._is3D());
        this._options.z.enable();

        const symbols = this._symbols();
        for (let s = 0; s < this._data.maxSymbols; s++) {
            symbols.push([get3DSymbol(s)]);
        }

        // switch all traces to 3D mode
        this._restyle({
            'marker.symbol': symbols,
            type: 'scatter3d',
        } as unknown as Data);

        for (const data of this._selected.values()) {
            data.toggleVisible(false);
        }
        this._updateMarkers();

        // Change the data that vary between 2D and 3D mode
        this._restyle(
            {
                // transparency messes with depth sorting in 3D mode, even with
                // line width set to 0 ¯\_(ツ)_/¯
                // https://github.com/plotly/plotly.js/issues/4111
                'marker.line.color': this._lineColors(),
                'marker.line.width': [1, 2],
                // size change from 2D to 3D
                'marker.size': this._sizes(),
                'marker.sizemode': 'area',
            } as Data,
            [0, 1]
        );

        this._relayout({
            // change colorbar length to accomodate for symbols legend
            'coloraxis.colorbar.len': this._colorbarLen(),
            // Carry over axis types
            'scene.xaxis.type': this._options.x.scale.value as Plotly.AxisType,
            'scene.yaxis.type': this._options.y.scale.value as Plotly.AxisType,
            'scene.zaxis.type': this._options.z.scale.value as Plotly.AxisType,
        } as unknown as Layout);
    }

    /** Switch current plot from 3D back to 2D */
    private _switch2D(): void {
        assert(!this._is3D());
        this._options.z.disable();

        const symbols = this._symbols();
        for (let sym = 0; sym < this._data.maxSymbols; sym++) {
            symbols.push([sym]);
        }

        // switch all traces to 2D mode
        this._restyle({
            'marker.symbol': symbols,
            type: 'scattergl',
        } as unknown as Data);

        // show selected environments markers
        for (const data of this._selected.values()) {
            data.toggleVisible(true);
        }

        this._restyle(
            {
                // transparency messes with depth sorting in 3D mode
                // https://github.com/plotly/plotly.js/issues/4111
                'marker.line.color': this._lineColors(),
                'marker.line.width': [1, 0],
                // size change from 2D to 3D
                'marker.size': this._sizes(),
            } as Data,
            [0, 1]
        );

        this._relayout({
            // change colorbar length to accomodate for symbols legend
            'coloraxis.colorbar.len': this._colorbarLen(),
            // Carry over axis types
            'xaxis.type': this._options.x.scale.value as Plotly.AxisType,
            'yaxis.type': this._options.y.scale.value as Plotly.AxisType,
        } as unknown as Layout);
    }

    /**
     * Function used as callback to update the axis ranges in settings after
     * the user changes zoom or range on the plot
     */
    private _afterplot(): void {
        // HACK: this is not public, so it might break
        const bounds = this._getRange();

        const xbound = bounds.get('x');
        const ybound = bounds.get('y');
        const zbound = bounds.get('z');

        assert(xbound !== undefined && ybound !== undefined);
        this._options.x.min.value = xbound[0];
        this._options.x.max.value = xbound[1];
        this._options.y.min.value = ybound[0];
        this._options.y.max.value = ybound[1];
        if (zbound !== undefined) {
            this._options.z.min.value = zbound[0];
            this._options.z.max.value = zbound[1];
        }

        if (!this._is3D()) {
            this._updateMarkers();
        }
    }

    /**
     * Update the position, color & size of markers within the data array
     */
    private _updateMarkers(data: MarkerData[] = Array.from(this._selected.values())): void {
        if (this._is3D()) {
            data.forEach((d) => d.toggleVisible(false));
            this._restyle(
                {
                    'marker.color': this._colors(1),
                    'marker.size': this._sizes(1),
                    'marker.symbol': this._symbols(1),
                    x: this._coordinates(this._options.x, 1),
                    y: this._coordinates(this._options.y, 1),
                    z: this._coordinates(this._options.z, 1),
                } as Data,
                1
            );
        } else {
            const allX = this._coordinates(this._options.x, 0) as number[][];
            const allY = this._coordinates(this._options.y, 0) as number[][];
            const plotWidth = this._plot.getBoundingClientRect().width;

            for (const datum of data) {
                const rawX = allX[0][datum.current];
                const rawY = allY[0][datum.current];
                if (this._insidePlot(rawX, rawY)) {
                    const x = plotWidth - this._computeRSCoord(rawX, 'x');
                    const y = this._computeRSCoord(rawY, 'y');
                    datum.update(x, y);
                } else {
                    datum.toggleVisible();
                }
            }
        }
    }

    // Get the current boundaries. Returns map of 'x' | 'y' | 'z' --> limits
    private _getRange(): Map<string, number[]> {
        const bounds = new Map<string, number[]>();
        let layout;
        if (this._is3D()) {
            layout = this._plot._fullLayout.scene;
            bounds.set('x', layout.xaxis.range as number[]);
            bounds.set('y', layout.yaxis.range as number[]);
            bounds.set('z', layout.zaxis.range as number[]);
        } else {
            layout = this._plot._fullLayout;
            bounds.set('x', layout.xaxis.range as number[]);
            bounds.set('y', layout.yaxis.range as number[]);
        }
        return bounds;
    }

    // Computes the real space coordinate of a value on the plot
    private _computeRSCoord(value: number, axisName: string): number {
        assert(axisName === 'x' || axisName === 'y');
        assert(!this._is3D());
        let axis;
        switch (axisName) {
            case 'x':
                axis = this._plot._fullLayout.xaxis;
                break;
            case 'y':
                axis = this._plot._fullLayout.yaxis;
                break;
        }
        return axis.l2p(value) + axis._offset;
    }

    // Checks if a point is in the visible plot for a *single* axis
    private _checkBounds(value: number, axis: string, buffer: number = 10): boolean {
        const bounds = this._getRange().get(axis);
        assert(bounds !== undefined);
        return value > bounds[0] - buffer && value < bounds[1] + buffer;
    }

    // Checks if a point is in the visible plot
    private _insidePlot(x: number, y: number, z?: number, buffer: number = 10): boolean {
        const check = this._checkBounds(x, 'x', buffer) && this._checkBounds(y, 'y', buffer);
        if (z !== undefined) {
            return this._checkBounds(z, 'z', buffer) && check;
        }
        return check;
    }
}
