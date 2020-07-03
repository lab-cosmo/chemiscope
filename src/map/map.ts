/**
 * @packageDocumentation
 * @module map
 */

import assert from 'assert';

import {Property} from '../dataset';
import {EnvironmentIndexer, HTMLSetting, Indexes, PositioningCallback, SettingGroup, SettingModificationOrigin} from '../utils';
import {foreachSetting, getByID, makeDraggable, sendWarning} from '../utils';

import Plotly from './plotly/plotly-scatter';
import {Config, Data, Layout, PlotlyScatterElement} from './plotly/plotly-scatter';

import {COLOR_MAPS} from './colorscales';
import {MapData, NumericProperty} from './data';

import BARS_SVG from '../static/bars.svg';
import HTML_SETTINGS from './settings.html';

const DEFAULT_LAYOUT = {
    // coloraxis is used for the markers
    coloraxis: {
        cmax: 0,
        cmin: 0,
        colorbar: {
            len: 1,
            thickness: 20,
            title: {
                text: '' as (undefined | string),
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
            title: '' as (undefined | string),
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
    ],
};

// in 3D mode, only strings are supported for 'marker.symbol', and only very few
// of them. See https://github.com/plotly/plotly.js/issues/4205 as the plotly
// issue tracking more symbols in 3D mode.
const POSSIBLE_SYMBOLS_IN_3D = ['circle', 'square', 'diamond', 'cross', 'x'];

function get3DSymbol(i: number): string {
    return POSSIBLE_SYMBOLS_IN_3D[i % POSSIBLE_SYMBOLS_IN_3D.length];
}

// get the max/min of an array. Math.min(...array) fails with very large arrays
function arrayMaxMin(values: number[]): {max: number, min: number} {
    let max = Number.NEGATIVE_INFINITY;
    let min = Number.POSITIVE_INFINITY;
    for (const value of values) {
        if (value > max) {
            max = value;
        }
        if (value < min) {
            min = value;
        }
    }
    assert(isFinite(min) && isFinite(max));
    return {max, min};
}

/** HTML element holding settings for a given axis (x, y, z, color) */
class AxisSetting {
    /// Which property should we use for this axis
    public property: HTMLSetting<'string'>;
    /// Which scale (linear/log) should we use
    public scale: HTMLSetting<'string'>;
    /// The minimal value for this axis
    public min: HTMLSetting<'number'>;
    /// The maximal value for this axis
    public max: HTMLSetting<'number'>;

    constructor(properties: string[], initial: string) {
        this.max = new HTMLSetting('number', 0);
        this.min = new HTMLSetting('number', 0);

        this.property = new HTMLSetting('string', initial);
        this.property.validate = (value) => {
            if (properties.includes(value) || (initial === '' && value === '')) {
                return;
            }
            throw Error(`invalid property '${value}' for axis`);
        };

        this.scale = new HTMLSetting('string', 'linear');
        this.scale.validate = (value) => {
            if (value === 'linear' || value === 'log') {
                return;
            }
            throw Error(`invalid value '${value}' for axis scale`);
        };
    }

    /** Disable auxiliary settings (min/max/scale) related to this axis */
    public disable() {
        this.max.disable();
        this.min.disable();
        this.scale.disable();
    }

    /** Enable auxiliary settings (min/max/scale) related to this axis */
    public enable() {
        this.max.enable();
        this.min.enable();
        this.scale.enable();
    }
}

/// interface to contain synchronized parameters for marker instances
interface MarkerData {
    /// color of the marker, synchronized with the StructureViewer where appropriate
    color: string;
    /// index of the currently displayed structure
    current: number;
    /// Marker indicating the position of the point in 2D mode
    /// Using an HTML element is much faster than using Plolty to restyle the full plot,
    /// especially with more than 100k points. In 3D mode, a separate trace is
    /// used instead.
    marker: HTMLElement;
}
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
    /** Callback fired when the plot is clicked and a new point is selected */
    public onselect: (indexes: Indexes) => void;
    /** Callback fired when a guid is removed from the map */
    public onremove: (removedGUID: string) => void;
    /** Callback fired when a new active marker is designated */
    public activate: (activeGUID: string) => void;

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
    private _active?: string;
    /// Map associating currently selected markers GUID to additional data
    private _selected: Map<string, MarkerData>;

    /// environment indexer
    private _indexer: EnvironmentIndexer;
    /// Store the HTML elements used for settings
    private _settings!: {
        x: AxisSetting;
        y: AxisSetting;
        z: AxisSetting;
        color: AxisSetting;
        palette: HTMLSetting<'string'>;
        symbol: HTMLSetting<'string'>;
        size: {
            property: HTMLSetting<'string'>;
            factor: HTMLSetting<'number'>;
        };
    };
    /// Button used to reset the range of color axis
    private _colorReset: HTMLButtonElement;
    /// The HTML element containing the settings modal
    private _settingsModal!: HTMLElement;

    /**
     * Create a new [[PropertiesMap]] inside the DOM element with the given HTML
     * `id`
     *
     * @param id         HTML id of the DOM element where the map should live
     * @param indexer    [[EnvironmentIndexer]] used to translate indexes from
     *                   environments index to structure/atom indexes
     * @param properties properties to be displayed
     */
    constructor(id: string,
                indexer: EnvironmentIndexer,
                properties: {[name: string]: Property},
                ) {
        this._indexer = indexer;
        this.onselect = () => {};
        this.onremove = () => {};
        this.activate = () => {};
        this._selected = new Map();

        this._root = getByID(id);

        if (this._root.style.position === '') {
            this._root.style.position = 'relative';
        }

        this._plot = document.createElement('div') as unknown as PlotlyScatterElement;
        this._plot.style.width = '100%';
        this._plot.style.height = '100%';
        this._root.appendChild(this._plot);

        this._data = new MapData(properties);

        this._insertSettingsHTML();
        this._colorReset = getByID<HTMLButtonElement>('chsp-color-reset');

        this._setupSettings();
        this._connectSettings();

        this._createPlot();

        // By default, position the modal for settings on top of the plot,
        // centered horizontally
        this.positionSettingsModal = (rect: DOMRect) => {
            const rootRect = this._root.getBoundingClientRect();
            return {
                left: rootRect.left + rootRect.width / 2 - rect.width / 2,
                top: rootRect.top + 20,
            };
        };
    }

    /** Change the environment corresponding to the `selectGUID` in
     * this._selected to the one with the given `indexes`
     * if selectedGUID is not defined, it will default to the active GUID, if available
     */
    public select(indexes: Indexes, selectedGUID?: string) {
        // Plotly.restyle fires the plotly_click event, so ensure we only run
        // the update once.
        // https://github.com/plotly/plotly.js/issues/1025

        if (selectedGUID === undefined) {
          if (this.active !== undefined) {
            selectedGUID = this.active;
          } else { return; }
        }

        /// Checks if marker exists on this map, if not adds it
        if (selectedGUID !== undefined && !this._selected.has(selectedGUID)) {
            this._addMarker(selectedGUID, indexes.environment);

            const newMarkerData = this._selected.get(selectedGUID);
            assert(newMarkerData !== undefined);
            this._updateSelectedMarker(selectedGUID, newMarkerData);
        }

        // sets this marker to active
        this.active = selectedGUID;

        const markerData = this._selected.get(selectedGUID);
        assert(markerData !== undefined);

        /// Sets the active marker on this map
        /// default of markerData.current is -1 so that currents of value 0 are
        /// still updated
        if (markerData.current === 0 || indexes.environment !== markerData.current ) {
            markerData.current = indexes.environment;
            this._updateSelectedMarker(selectedGUID, markerData);
        }

    }

    /**
     * Remove all HTML added by this [[PropertiesMap]] in the current document
     */
    public remove(): void {
        this._root.innerHTML = '';
        this._settingsModal.remove();
    }

    public get active(): string {
      if (this._selected.size > 0) {
        if (this._active === undefined) {
            this.active = this._selected.keys().next().value;
            return this._selected.keys().next().value;
        } else {
            return this._active; }
      } else {
            return ''; }
    }
    /**
     * Function to set the active marker for communicating with the structure
     * viewer
     *
     * @param  activeGUID the GUID of the new active viewer
     */
    public set active(activeGUID: string) {
        if (this._active === activeGUID) {
            // guard against infinite recursion
            // `Structure set active => StructureViewer.show => Map set active
            // => PropertiesMap.select => Structure set active`
            return;
        }

        if (!this._selected.has(activeGUID)) {
          this._addMarker(activeGUID);
        }

        if (!this._is3D()) {
            // if we are in 2D we need to update the visuals
            if (this._active !== undefined && this._selected.has(this._active)) {
                // this has been left as document.getElementById for when we are
                // setting a new active because the old one has been deleted.
                const oldButton = document.getElementById(`chsp-selected-${this._active}`);
                if (oldButton !== null) {
                    oldButton.classList.toggle('chsp-active-structure', false);
                }
            }
            const newButton = getByID(`chsp-selected-${activeGUID}`);
            this._active = activeGUID;
            newButton.classList.toggle('chsp-active-structure', true);

        } else {
            this._active = activeGUID;
            const factor = this._settings.size.factor.value;
            this._restyle({'marker.size': this._sizes(factor, 1)} as Data, 1);
        }

        const markerData = this._selected.get(this._active);
        assert(markerData !== undefined);

        this.activate(this._active);
    }

    /**
     * Removes a marker from the map.
     *
     * @param  removedGUID     unique string identifier of the marker to remove
     */
    public removeMarker(removedGUID: string): void {
        const marker = document.getElementById(`chsp-selected-${removedGUID}`);
        if (marker !== null) {
            if (marker.parentNode !== null) {
                marker.parentNode.removeChild(marker);
            }

            this._selected.delete(removedGUID);

            if (this._active === removedGUID) {
                assert(this._selected.size > 0);
                this.active = this._selected.keys().next().value;
            }
            this.onremove(removedGUID);
        }
    }

    /**
     * Function to hide the marker when switching from 2D to 3D
     */
    private _hideMarker(hiddenGUID: string): void {
      const marker = document.getElementById(`chsp-selected-${hiddenGUID}`);
      if (marker !== null) {
          if (marker.parentNode !== null) {
              marker.parentNode.removeChild(marker);
          }
      }
    }

    /** Forward to Ploty.restyle */
    private _restyle(data: Partial<Data>, traces?: number | number[]) {
        Plotly.restyle(this._plot, data, traces).catch((e) => setTimeout(() => { throw e; }));
    }

    /** Forward to Ploty.relayout */
    private _relayout(layout: Partial<Layout>) {
        Plotly.relayout(this._plot, layout).catch((e) => setTimeout(() => { throw e; }));
    }

    /** Create the settings modal by adding HTML to the page */
    private _insertSettingsHTML() {
        // use HTML5 template to generate a DOM object from an HTML string
        const template = document.createElement('template');
        template.innerHTML = `<button
            class="btn btn-light btn-sm chsp-viewer-button"
            data-target='#chsp-settings'
            data-toggle="modal"
            style="top: 4px; left: 5px; opacity: 1;">
                <div>${BARS_SVG}</div>
            </button>`;
        const openSettings = template.content.firstChild!;
        this._root.append(openSettings);

        // TODO: set unique HTML id in the settings to allow multiple map in
        // the same page
        template.innerHTML = HTML_SETTINGS;
        this._settingsModal = template.content.firstChild! as HTMLElement;
        document.body.appendChild(this._settingsModal);

        const modalDialog = this._settingsModal.childNodes[1]! as HTMLElement;
        if (!modalDialog.classList.contains('modal-dialog')) {
            throw Error('internal error: missing modal-dialog class');
        }
        // make the settings modal draggable
        makeDraggable(modalDialog, '.modal-header');

        // Position modal near the actual viewer
        openSettings.addEventListener('click', () => {
            // only set style once, on first open, and keep previous position
            // on next open to keep the 'draged-to' position
            if (modalDialog.getAttribute('data-initial-modal-positions-set') === null) {
                modalDialog.setAttribute('data-initial-modal-positions-set', 'true');

                // display: block to ensure modalDialog.offsetWidth is non-zero
                (modalDialog.parentNode as HTMLElement).style.display = 'block';

                const {top, left} = this.positionSettingsModal(modalDialog.getBoundingClientRect());

                // set width first, since setting position can influence it
                modalDialog.style.width = `${modalDialog.offsetWidth}px`;
                // unset margins when using position: fixed
                modalDialog.style.margin = '0';
                modalDialog.style.position = 'fixed';
                modalDialog.style.top = `${top}px`;
                modalDialog.style.left = `${left}px`;
            }
        });
    }

    /** Add all the required callback to the settings */
    private _connectSettings() {
        // ======= x axis settings
        this._settings.x.property.onchange = () => {
            const values = this._xValues();
            this._restyle({ x: values }, [0, 1]);
            this._relayout({
                'scene.xaxis.title': this._settings.x.property.value,
                'xaxis.title': this._settings.x.property.value,
            } as unknown as Layout);
        };

        this._settings.x.scale.onchange = () => {
            if (this._is3D()) {
                this._relayout({
                    'scene.xaxis.type': this._settings.x.scale.value,
                } as unknown as Layout);
            } else {
                this._relayout({ 'xaxis.type': this._settings.x.scale.value as Plotly.AxisType });
            }
        };

        // function creating a function to be used as onchange callback
        // for <axis>.min and <axis>.max
        const rangeChange = (name: string, axis: AxisSetting) => {
            return (_: number, origin: SettingModificationOrigin) => {
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

        this._settings.x.min.onchange = rangeChange('xaxis', this._settings.x);
        this._settings.x.max.onchange = rangeChange('xaxis', this._settings.x);

        // ======= y axis settings
        this._settings.y.property.onchange = () => {
            const values = this._yValues();
            this._restyle({ y:  values}, [0, 1]);
            this._relayout({
                'scene.yaxis.title': this._settings.y.property.value,
                'yaxis.title': this._settings.y.property.value,
            } as unknown as Layout);
        };

        this._settings.y.scale.onchange = () => {
            if (this._is3D()) {
                this._relayout({
                    'scene.yaxis.type': this._settings.y.scale.value,
                } as unknown as Layout);
            } else {
                this._relayout({ 'yaxis.type': this._settings.y.scale.value as Plotly.AxisType });
            }
        };

        this._settings.y.min.onchange = rangeChange('yaxis', this._settings.y);
        this._settings.y.max.onchange = rangeChange('yaxis', this._settings.y);

        // ======= z axis settings
        this._settings.z.property.onchange = () => {
            const was3D = (this._plot as any)._fullData[0].type === 'scatter3d';
            if (this._settings.z.property.value === '') {
                if (was3D) {
                    this._switch2D();
                }
            } else {
                if (!was3D) {
                    this._switch3D();
                }
            }

            const values = this._zValues();
            this._restyle({ z: values } as Data, [0, 1]);
            this._relayout({
                'scene.zaxis.title': this._settings.z.property.value,
            } as unknown as Layout);
        };

        this._settings.z.scale.onchange = () => {
            this._relayout({
                'scene.zaxis.type': this._settings.z.scale.value,
            } as unknown as Layout);
        };

        this._settings.z.min.onchange = rangeChange('zaxis', this._settings.z);
        this._settings.z.max.onchange = rangeChange('zaxis', this._settings.z);

        // ======= color axis settings
        this._settings.color.property.onchange = () => {
            if (this._settings.color.property.value !== '') {
                this._settings.color.enable();
                this._colorReset.disabled = false;

                const values = this._colors(0)[0] as number[];
                const {min, max} = arrayMaxMin(values);

                this._settings.color.min.value = min;
                this._settings.color.max.value = max;

                this._relayout({
                    'coloraxis.colorbar.title.text': this._settings.color.property.value,
                    'coloraxis.showscale': true,
                } as unknown as Layout);

            } else {
                this._settings.color.disable();
                this._colorReset.disabled = true;

                this._settings.color.min.value = 0;
                this._settings.color.max.value = 0;

                this._relayout({
                    'coloraxis.colorbar.title.text': undefined,
                    'coloraxis.showscale': false,
                } as unknown as Layout);
            }

            this._restyle({
                'hovertemplate': this._hovertemplate(),
                'marker.color': this._colors(0),
            } as Data, 0);
        };

        const colorRangeChange = () => {
            const min = this._settings.color.min.value;
            const max = this._settings.color.max.value;
            this._relayout({
                'coloraxis.cmax': max,
                'coloraxis.cmin': min,
                // looks like changing only 'coloraxis.cmax'/'coloraxis.cmin' do
                // not update the color of the points (although it does change
                // the colorbar). Asking for an update of 'coloraxis.colorscale'
                // seems to do the trick. This is possiblely a Ploty bug, we
                // would need to investiguate a bit more.
                'coloraxis.colorscale': this._colorScale(),
            } as unknown as Layout);
        };
        this._settings.color.min.onchange = colorRangeChange;
        this._settings.color.max.onchange = colorRangeChange;

        this._colorReset.onclick = () => {
            const values = this._colors(0)[0] as number[];
            const {min, max} = arrayMaxMin(values);
            this._settings.color.min.value = min;
            this._settings.color.max.value = max;
            this._relayout({
                'coloraxis.cmax': max,
                'coloraxis.cmin': min,
                // same as above regarding update of the points color
                'coloraxis.colorscale': this._colorScale(),
            } as unknown as Layout);
        };

        // ======= color palette
        this._settings.palette.onchange = () => {
            this._relayout({
                'coloraxis.colorscale': this._colorScale(),
            } as unknown as Layout);
        };

        // ======= markers symbols
        this._settings.symbol.onchange = () => {
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
        this._settings.size.property.onchange = () => {
            const factor = this._settings.size.factor.value;
            this._restyle({ 'marker.size': this._sizes(factor, 0) } as Data, 0);
        };

        this._settings.size.factor.onchange = () => {
            const factor = this._settings.size.factor.value;
            this._restyle({ 'marker.size': this._sizes(factor, 0) } as Data, 0);
        };
    }

    /**
     * Fill possible values for the settings, depending on the properties in
     * the dataset
     */
    private _setupSettings() {
        const properties = Object.keys(this._properties());
        if (properties.length < 2) {
            throw Error('we need at least two properties to plot in the map');
        }

        if (this._settings !== undefined) {
            // when changing dataset, remove all previous event listeners
            foreachSetting(this._settings as unknown as SettingGroup, (setting) => {
                setting.unbindAll();
            });
        }

        this._settings = {
            color: new AxisSetting(properties, ''),
            x: new AxisSetting(properties, properties[0]),
            y: new AxisSetting(properties, properties[1]),
            z: new AxisSetting(properties, ''),

            palette: new HTMLSetting('string', 'inferno'),
            size: {
                factor: new HTMLSetting('number', 50),
                property: new HTMLSetting('string', ''),
            },
            symbol: new HTMLSetting('string', ''),
        };

        const validate = (value: string) => {
            if (properties.includes(value) || value === '') {
                return;
            }
            throw Error(`invalid property name '${value}'`);
        };
        this._settings.size.property.validate = validate;
        this._settings.symbol.validate = validate;

        // ============== Setup the map options ==============
        // ======= data used as x values
        const selectXProperty = getByID<HTMLSelectElement>('chsp-x');
        selectXProperty.options.length = 0;
        for (const key of properties) {
            selectXProperty.options.add(new Option(key, key));
        }
        this._settings.x.property.bind(selectXProperty, 'value');
        this._settings.x.min.bind('chsp-x-min', 'value');
        this._settings.x.max.bind('chsp-x-max', 'value');
        this._settings.x.scale.bind('chsp-x-scale', 'value');

        // ======= data used as y values
        const selectYProperty = getByID<HTMLSelectElement>('chsp-y');
        selectYProperty.options.length = 0;
        for (const key of properties) {
            selectYProperty.options.add(new Option(key, key));
        }
        this._settings.y.property.bind(selectYProperty, 'value');
        this._settings.y.min.bind('chsp-y-min', 'value');
        this._settings.y.max.bind('chsp-y-max', 'value');
        this._settings.y.scale.bind('chsp-y-scale', 'value');

        // ======= data used as z values
        const selectZProperty = getByID<HTMLSelectElement>('chsp-z');
        // first option is 'none'
        selectZProperty.options.length = 0;
        selectZProperty.options.add(new Option('none', ''));
        for (const key of properties) {
            selectZProperty.options.add(new Option(key, key));
        }
        this._settings.z.property.bind(selectZProperty, 'value');
        this._settings.z.min.bind('chsp-z-min', 'value');
        this._settings.z.max.bind('chsp-z-max', 'value');
        this._settings.z.scale.bind('chsp-z-scale', 'value');

        // ======= data used as color values
        const selectColorProperty = getByID<HTMLSelectElement>('chsp-color');
        // first option is 'none'
        selectColorProperty.options.length = 0;
        selectColorProperty.options.add(new Option('none', ''));
        for (const key of properties) {
            selectColorProperty.options.add(new Option(key, key));
        }
        this._settings.color.property.bind(selectColorProperty, 'value');
        this._settings.color.min.bind('chsp-color-min', 'value');
        this._settings.color.max.bind('chsp-color-max', 'value');

        if (properties.length >= 3) {
            this._settings.color.property.value = properties[2];
            this._settings.color.enable();
            this._colorReset.disabled = false;

            const values = this._colors(0)[0] as number[];
            const {min, max} = arrayMaxMin(values);

            this._settings.color.min.value = min;
            this._settings.color.max.value = max;
        } else {
            this._settings.color.property.value = '';
            this._settings.color.min.disable();
            this._colorReset.disabled = true;

            this._settings.color.min.value = 0;
            this._settings.color.max.value = 0;
        }

        // ======= color palette
        const selectPalette = getByID<HTMLSelectElement>('chsp-palette');
        selectPalette.length = 0;
        for (const key in COLOR_MAPS) {
            selectPalette.options.add(new Option(key, key));
        }
        this._settings.palette.bind(selectPalette, 'value');

        // ======= marker symbols
        const selectSymbolProperty = getByID<HTMLSelectElement>('chsp-symbol');
        // first option is 'default'
        selectSymbolProperty.options.length = 0;
        selectSymbolProperty.options.add(new Option('default', ''));
        for (const key of properties) {
            if (this._property(key).string !== undefined) {
                selectSymbolProperty.options.add(new Option(key, key));
            }
        }
        this._settings.symbol.bind(selectSymbolProperty, 'value');

        // ======= marker size
        const selectSizeProperty = getByID<HTMLSelectElement>('chsp-size');
        // first option is 'default'
        selectSizeProperty.options.length = 0;
        selectSizeProperty.options.add(new Option('default', ''));
        for (const key of properties) {
            selectSizeProperty.options.add(new Option(key, key));
        }
        this._settings.size.property.bind(selectSizeProperty, 'value');
        this._settings.size.factor.bind('chsp-size-factor', 'value');
    }

    /** Actually create the Plotly plot */
    private _createPlot() {
        this._plot.innerHTML = '';

        const colors = this._colors();
        const lineColors = this._lineColors();
        // default value for the size factor is 50
        const sizes = this._sizes(50);
        const symbols = this._symbols();

        const x = this._xValues();
        const y = this._yValues();
        const z = this._zValues();

        // The main trace, containing default data
        const main = {
            name: '',
            type: 'scattergl',

            x: x[0],
            y: y[0],
            z: z[0],

            hovertemplate: this._hovertemplate(),
            marker: {
                color: colors[0],
                coloraxis: 'coloraxis',
                line: {
                    color: lineColors[0],
                    width: 1,
                },
                // prevent plolty from messing with opacity when doing bubble
                // style charts (different sizes for each point)
                opacity: 1,
                size: sizes[0],
                sizemode: 'area',
                symbol: symbols[0],
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
            type: 'scattergl',

            x: [NaN],
            y: [NaN],
            z: [NaN],

            hoverinfo: 'none',
            marker: {
                color: colors[1],
                line: {
                    color: lineColors[1],
                    width: 2,
                },
                opacity: 1,
                size: sizes[1],
                sizemode: 'area',
            },
            mode: 'markers',
            showlegend: false,
        };

        const traces = [main as Data, selected as Data];

        // add empty traces to be able to display the symbols legend
        // one trace for each possible symbol
        for (let s = 0; s < this._data.maxSymbols; s++) {
            const data = {
                name: '',
                type: 'scattergl',

                x: [NaN],
                y: [NaN],
                z: [NaN],

                marker: {
                    color: 'black',
                    size: 10,
                    symbol: s,
                },
                mode: 'markers',
                showlegend: false,
            };
            traces.push(data as Data);
        }

        // make a copy of the default layout
        const layout = JSON.parse(JSON.stringify(DEFAULT_LAYOUT));
        // and set values speific to the displayed dataset
        layout.xaxis.title = this._settings.x.property.value;
        layout.yaxis.title = this._settings.y.property.value;
        layout.xaxis.type = this._settings.x.scale.value;
        layout.yaxis.type = this._settings.y.scale.value;
        layout.scene.xaxis.title = this._settings.x.property.value;
        layout.scene.yaxis.title = this._settings.y.property.value;
        layout.scene.zaxis.title = this._settings.z.property.value;
        layout.coloraxis.colorscale = this._colorScale();
        layout.coloraxis.cmin = this._settings.color.min.value;
        layout.coloraxis.cmax = this._settings.color.max.value;
        layout.coloraxis.colorbar.title.text = this._settings.color.property.value;

        // Create an empty plot and fill it below
        Plotly.newPlot(this._plot, traces, layout as Partial<Layout>, DEFAULT_CONFIG as Config)
            .catch((e) => setTimeout(() => { throw e; }));

        this._plot.on('plotly_click', (event: Plotly.PlotMouseEvent) => {
            // don't update selected env on double click, since it is lined to
            // 'reset zoom level' in 2D mode.
            // `event.event` is only set in 2D mode

            if (event.event && event.event.detail === 2) {
                return;
            }

            let data;
            // if someone has clicked on a selection marker, set to active
            if (event.points[0].data.name === 'selected' ) {
                let i = 0;
                for (const [guid, markerData] of this._selected.entries()) {
                    if (i === event.points[0].pointNumber) {
                        this.active = guid;
                        data = markerData;
                        break;
                    }
                    i++;
                }
            }

            const indexes = (data !== undefined) ?
                this._indexer.from_environment(data.current) :
                // default if someone has clicked generally on the map or on a
                // place the selected marker doesn't recognize
                this._indexer.from_environment(event.points[0].pointNumber);

            this.select(indexes, this._active);
            this.onselect(indexes);
        });
        this._plot.on('plotly_afterplot', () => this._afterplot());

        this._updateAllMarkers();
    }

    /** Get the currently available properties: either `'atom'` or `'structure'` properties */
    private _properties(): {[name: string]: NumericProperty} {
        return this._data[this._indexer.mode];
    }

    /** Get the property with the given name */
    private _property(name: string): NumericProperty {
        const result = this._data[this._indexer.mode][name];
        if (result === undefined) {
            throw Error(`unknown property '${name}' requested in map`);
        }
        return result;
    }

    /** Get the plotly hovertemplate depending on `this._current.color` */
    private _hovertemplate(): string {
        if (this._hasColors()) {
            return this._settings.color.property.value + ': %{marker.color:.2f}<extra></extra>';
        } else {
            return '%{x:.2f}, %{y:.2f}<extra></extra>';
        }
    }

    /**
     * Get the values to use for the x axis with the given plotly `trace`,
     * or all of them if `trace === undefined`
     */
    private _xValues(trace?: number): number[][] {
        const values = this._property(this._settings.x.property.value).values;
        const selected = [];
        for (const marker of this._selected.values()) {
            if (this._is3D()) {
                selected.push(values[marker.current]);
            } else {
                selected.push(NaN);
            }
        }
        return this._selectTrace(values, selected, trace);
    }

    /**
     * Get the values to use for the y axis with the given plotly `trace`,
     * or all of them if `trace === undefined`
     */
    private _yValues(trace?: number): number[][] {
        const values = this._property(this._settings.y.property.value).values;
        const selected = [];
        for (const marker of this._selected.values()) {
            if (this._is3D()) {
                selected.push(values[marker.current]);
            } else {
                selected.push(NaN);
            }
        }

        return this._selectTrace(values, selected, trace);
    }

    /**
     * Get the values to use for the z axis with the given plotly `trace`,
     * or all of them if `trace === undefined`
     */
    private _zValues(trace?: number): Array<undefined | number[]> {
        if (!this._is3D()) {
            return this._selectTrace(undefined, undefined, trace);
        }

        const values = this._property(this._settings.z.property.value).values;
        const selected = [];
        for (const marker of this._selected.values()) {
            selected.push(values[marker.current]);
        }
        return this._selectTrace(values, selected, trace);
    }

    /**
     * Get the color values to use with the given plotly `trace`, or all of
     * them if `trace === undefined`
     */
    private _colors(trace?: number): Array<string | string[] | number | number[]> {
        let values;
        if (this._hasColors()) {
            values = this._property(this._settings.color.property.value).values;
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

    /** Get the colorscale to use for markers in the main plotly trace */
    private _colorScale(): Plotly.ColorScale {
        return COLOR_MAPS[this._settings.palette.value];
    }

    /**
     * Get the values to use as marker size with the given plotly `trace`, or
     * all of them if `trace === undefined`.
     *
     * The size scaling parameter should be given in `sizeSliderValue`.
     */
    private _sizes(sizeSliderValue: number, trace?: number): Array<number | number[]> {
        // Transform the linear value from the slider into a logarithmic scale
        const logSlider = (value: number) => {
            const min_slider = 1;
            const max_slider = 100;

            // go from 1/6th of the size to 6 time the size
            const min_value = Math.log(1.0 / 6.0);
            const max_value = Math.log(6.0);

            const tmp = (max_value - min_value) / (max_slider - min_slider);
            return Math.exp(min_value + tmp * (value - min_slider));
        };

        const userFactor = logSlider(sizeSliderValue);

        let values;
        if (this._settings.size.property.value !== '') {
            const sizes = this._property(this._settings.size.property.value).values;
            const {min, max} = arrayMaxMin(sizes);
            const defaultSize = this._is3D() ? 2000 : 150;
            values = sizes.map((v: number) => {
                // normalize between 0 and 1, then scale by the user provided value
                const scaled = userFactor * (v + 0.05 - min) / (max - min);
                // since we are using scalemode: 'area', square the scaled value
                return defaultSize * scaled * scaled;
            });
        } else {
            // we need to use an array instead of a single value because of
            // https://github.com/plotly/plotly.js/issues/2735
            values = new Array(this._indexer.environmentsCount());
            if (this._is3D()) {
                values.fill(500 * userFactor);
            } else {
                values.fill(50 * userFactor);
            }
        }

        const selected = [];
        if (this._is3D()) {
            for (const guid of this._selected.keys()) {
                if (guid === this._active) {
                    selected.push(4000);
                } else {
                    selected.push(2000);
                }
            }
        }
        return this._selectTrace<number | number[]>(values, selected, trace);
    }

    /**
     * Get the values to use as marker symbol with the given plotly `trace`, or
     * all of them if `trace === undefined`.
     */
    private _symbols(trace?: number): Array<number | number[] | string | string[]> {
        if (this._settings.symbol.value === '') {
            // default to 0 (i.e. circles)
            return this._selectTrace<string | string[]>('circle', 'circle', trace);
        }

        const values = this._property(this._settings.symbol.value).values;
        const markerData = this._selected.get(this.active);
        assert(markerData !== undefined);

        if (this._is3D()) {
            // If we need more symbols than available, we'll send a warning
            // and repeat existing ones
            if (this._symbolsCount() > POSSIBLE_SYMBOLS_IN_3D.length) {
              sendWarning(`${this._symbolsCount()} symbols are required, but we only have ${POSSIBLE_SYMBOLS_IN_3D.length}. Some symbols will be repeated`);
            }
            const symbols = values.map(get3DSymbol);
            const symbol = symbols[markerData.current];
            return this._selectTrace<string | string[]>(symbols, symbol, trace);
        } else {
            // in 2D mode, use automatic assignment of symbols from numeric
            // values
            const value = values[markerData.current];
            return this._selectTrace<number | number[]>(values, value, trace);
        }
    }

    /** Should we show the legend for the various symbols used? */
    private _showlegend(): boolean[] {
        const result = [false, false];

        if (this._settings.symbol.value !== '') {
            for (let i = 0; i < this._symbolsCount()!; i++) {
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

        if (this._settings.symbol.value !== '') {
            const names = this._property(this._settings.symbol.value).string!.strings();
            for (const name of names) {
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

    /** How many different symbols are being displayed */
    private _symbolsCount(): number {
        if (this._settings.symbol.value !== '') {
            return this._property(this._settings.symbol.value).string!.strings().length;
        } else {
            return 0;
        }
    }

    /** Get the length of the colorbar to accomodate for the legend */
    private _colorbarLen(): number {
        /// Heigh of a legend item in plot unit
        const LEGEND_ITEM_HEIGH = 0.045;
        return 1 - LEGEND_ITEM_HEIGH * this._symbolsCount();
    }

    /** Does the current plot use color values? */
    private _hasColors(): boolean {
        return this._settings.color.property.value !== '';
    }

    /** Is the the current plot a 3D plot? */
    private _is3D(): boolean {
        return this._settings !== undefined && this._settings.z.property.value !== '';
    }

    /** Switch current plot from 2D to 3D */
    private _switch3D() {
        assert(this._is3D());
        this._settings.z.enable();

        const symbols = this._symbols();
        for (let s = 0; s < this._data.maxSymbols; s++) {
            symbols.push([get3DSymbol(s)]);
        }

        // switch all traces to 3D mode
        this._restyle({
            'marker.symbol': symbols,
            'type': 'scatter3d',
        } as unknown as Data);

        const cachedActive = this.active;
        for (const [guid, markerData] of this._selected.entries()) {
            this._hideMarker(guid);
            this._updateSelectedMarker(guid, markerData);
        }
        this.active = cachedActive;

        // Change the data that vary between 2D and 3D mode
        const factor = this._settings.size.factor.value;
        this._restyle({
            // transparency messes with depth sorting in 3D mode, even with
            // line width set to 0 ¯\_(ツ)_/¯
            // https://github.com/plotly/plotly.js/issues/4111
            'marker.line.color': this._lineColors(),
            'marker.line.width': [0, 2],
            // size change from 2D to 3D
            'marker.size': this._sizes(factor),
            'marker.sizemode': 'area',
        } as Data, [0, 1]);

        this._relayout({
            // change colorbar length to accomodate for symbols legend
            'coloraxis.colorbar.len': this._colorbarLen(),
            // Carry over axis types
            'scene.xaxis.type': this._settings.x.scale.value as Plotly.AxisType,
            'scene.yaxis.type': this._settings.y.scale.value as Plotly.AxisType,
            'scene.zaxis.type': this._settings.z.scale.value as Plotly.AxisType,
        } as unknown as Layout);
    }

    /** Switch current plot from 3D back to 2D */
    private _switch2D() {
        assert(!this._is3D());
        this._settings.z.disable();

        const symbols = this._symbols();
        for (let sym = 0; sym < this._data.maxSymbols; sym++) {
            symbols.push([sym]);
        }

        // switch all traces to 2D mode
        this._restyle({
            'marker.symbol': symbols,
            'type': 'scattergl',
        } as unknown as Data);

        const cachedActive = this.active;

        for (const [guid, marker] of this._selected.entries()) {
            this._addMarker(guid, marker.current);
        }

        this.active = cachedActive;

        // Change the data that vary between 2D and 3D mode
        const factor = this._settings.size.factor.value;
        this._restyle({
            // transparency messes with depth sorting in 3D mode
            // https://github.com/plotly/plotly.js/issues/4111
            'marker.line.color': this._lineColors(),
            'marker.line.width': [1, 2],
            // size change from 2D to 3D
            'marker.size': this._sizes(factor),
        } as Data, [0, 1]);

        this._relayout({
            // change colorbar length to accomodate for symbols legend
            'coloraxis.colorbar.len': this._colorbarLen(),
            // Carry over axis types
            'xaxis.type': this._settings.x.scale.value as Plotly.AxisType,
            'yaxis.type': this._settings.y.scale.value as Plotly.AxisType,
        } as unknown as Layout);
    }

    /**
     * Function used as callback to update the axis ranges in settings after
     * the user changes zoom or range on the plot
     */
    private _afterplot() {
        // HACK: this is not public, so it might break
        const layout = this._plot._fullLayout;
        if (this._is3D()) {
            this._settings.x.min.value = layout.scene.xaxis.range[0];
            this._settings.x.max.value = layout.scene.xaxis.range[1];

            this._settings.y.min.value = layout.scene.yaxis.range[0];
            this._settings.y.max.value = layout.scene.yaxis.range[1];

            this._settings.z.min.value = layout.scene.zaxis.range[0];
            this._settings.z.max.value = layout.scene.zaxis.range[1];
        } else {
            this._settings.x.min.value = layout.xaxis.range[0];
            this._settings.x.max.value = layout.xaxis.range[1];

            this._settings.y.min.value = layout.yaxis.range[0];
            this._settings.y.max.value = layout.yaxis.range[1];

            this._updateAllMarkers();
        }
    }

    /**
     * Update the position of all the markers indicating selected points in the
     * map.
     *
     * In 3D mode, the markers uses the second Plotly trace.
     * In 2D mode, these markers are HTML div styled as colored circles that
     * we manually move around, saving a call to `restyle`.
     *
     * @param selectedGUID unique string identifier of the marker to update
     */
    private _updateSelectedMarker(selectedGUID: string, selectedMarkerData: MarkerData): void {
          const selected = selectedMarkerData.current;
          const marker = selectedMarkerData.marker;

          if (this._is3D()) {
              marker.style.display = 'none';
              let symbols;
              if (this._settings.symbol.value !== '') {
                  const values = this._property(this._settings.symbol.value).values;
                  symbols = get3DSymbol(values[selected]);
              } else {
                  symbols = get3DSymbol(0);
              }

              const factor = this._settings.size.factor.value;
              this._restyle({
                  'marker.color': this._colors(1),
                  'marker.size': this._sizes(factor, 1),
                  'marker.symbol': symbols,
                  'x': this._xValues(1),
                  'y': this._yValues(1),
                  'z': this._zValues(1),
              } as Data, 1);
          } else {
              const xaxis = this._plot._fullLayout.xaxis;
              const yaxis = this._plot._fullLayout.yaxis;

              const computeX = (data: number) => xaxis.l2p(data) + xaxis._offset;
              const computeY = (data: number) => yaxis.l2p(data) + yaxis._offset;

              const x = computeX(this._xValues(0)[0][selected]);
              const y = computeY(this._yValues(0)[0][selected]);

              // hide the point if it is outside the plot, allow for up to 10px
              // overflow (for points just on the border)
              const minX = computeX(xaxis.range[0]) - 10;
              const maxX = computeX(xaxis.range[1]) + 10;
              const minY = computeY(yaxis.range[1]) - 10;
              const maxY = computeY(yaxis.range[0]) + 10;
              if (!isFinite(x) || !isFinite(y) || x < minX || x > maxX || y < minY || y > maxY) {
                  marker.style.display = 'none';
              } else {
                  marker.style.display = 'block';
              }

              marker.style.top = `${y}px`;
              const plotWidth = this._plot.getBoundingClientRect().width;
              marker.style.right = `${plotWidth - x}px`;
          }
          const indexes = this._indexer.from_environment(selected);
          this.onselect(indexes);
    }

    private _updateAllMarkers() {
        // check if there are any markers to update
        if (this._selected.size > 0) {
          // stores info on current active point
          const active = this.active;
          const activeMarkerData = this._selected.get(this.active);
          assert(activeMarkerData !== undefined);
          const indexes = this._indexer.from_environment(activeMarkerData.current);

          for (const [guid, markerData] of this._selected.entries()) {
            this._updateSelectedMarker(guid, markerData);
          }

          // HACK: restores the point that was active before the update
          this.select(indexes, active);
        }
    }
    /**
     * Function to add a marker with the given GUID string and indices to the map.
     *
     * @param  addedGUID unique string identifier of the marker to add
     * @param  index     numeric index of the structure (with respect to dataset) to show
     */
    private _addMarker(addedGUID: string, index: number= 0): void {
        if (!this._selected.has(addedGUID)) {
            const activeButton = getByID(`chsp-activate-${addedGUID}`);
            const marker = document.createElement('div');
            marker.classList.add('chsp-structure-marker');
            if (addedGUID === this._active) {
                marker.classList.toggle('chsp-active-structure', true);
            }
            marker.id = `chsp-selected-${addedGUID}`;
            marker.onclick = () => this.active = addedGUID;

            const color = activeButton.style.backgroundColor;
            marker.style.backgroundColor = color;
            const newMarkerData = {
                                      color : color,
                                      current: index,
                                      marker: marker,
                                  };
            this._selected.set(addedGUID, newMarkerData);
        }

        const markerData = this._selected.get(addedGUID);
        assert(markerData !== undefined);

        this._root.appendChild(markerData.marker);

        if (this._is3D()) {
          markerData.marker.style.display = 'none';
        }
    }
}
