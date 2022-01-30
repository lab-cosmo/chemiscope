import { DOMWidgetView } from '@jupyter-widgets/base';
import { addWarningHandler, generateGUID, getByID } from '../../../src/utils';

// Import the CSS
import './widget.css';
import './chemiscope-bootstrap.less';
import 'bootstrap/dist/js/bootstrap.min.js';

import {
    DefaultVisualizer,
    MapVisualizer,
    StructureVisualizer,
    Settings,
} from '../../../src/index';
import { Dataset } from '../../../src/dataset';
import { SavedSettings } from '../../../src/options';

/**
 * The [[ChemiscopeView]] class renders the Chemiscope App as a widget in the
 * Jupyter Notebook output window when instantiated from the Chemiscope Python
 * package.
 */
export class ChemiscopeView extends DOMWidgetView {
    private visualizer?: DefaultVisualizer;
    private guid!: string;

    public render(): void {
        this.guid = `chsp-${generateGUID()}`;

        // this function works by first rendering the widget inside `this.el`,
        // and then inserting this.el inside the HTML document.
        const element = this.el;

        // HACK: resize Plotly when inserted. It currently render itself at full
        // width and then needs to be resized. For more on this, see
        // https://github.com/lab-cosmo/chemiscope/pull/181#discussion_r693005307
        element.addEventListener('DOMNodeInserted', () => {
            window.dispatchEvent(new Event('resize'));
        });

        addWarningHandler((message) => {
            const display = getByID(`${this.guid}-warning-display`, element);
            display.style.display = 'block';
            display.getElementsByTagName('p')[0].innerText = message;
        });

        element.innerHTML = `
        <div class="chemiscope-bootstrap">
            <div class="alert alert-warning" role="alert" id="${this.guid}-warning-display" style="display: none; font-size: 1.5em;">
                <button type="button" class="close" onclick="document.getElementById('${this.guid}-warning-display').style.display = 'none';">
                    <span aria-hidden="true">&times;</span>
                </button>
                <p></p>
            </div>
            <div class="alert alert-danger" role="alert" id="${this.guid}-error-display" style="display: none; font-size: 1.5em;">
                <button type="button" class="close" onclick="document.getElementById('${this.guid}-error-display').style.display = 'none';">
                    <span aria-hidden="true">&times;</span>
                </button>
                <p></p>
            </div>

            <div class="chemiscope-widget-two-col">
                <div class="chemiscope-meta-and-map">
                    <div class="chemiscope-meta" id="${this.guid}-chemiscope-meta"></div>
                    <div class="chemiscope-map" id="${this.guid}-chemiscope-map"></div>
                </div>
                <div class="chemiscope-structure-and-info">
                    <div class="chemiscope-structure" id="${this.guid}-chemiscope-structure"></div>
                    <div class="chemiscope-info" id="${this.guid}-chemiscope-info"></div>
                </div>
                </div>
            </div>
        </div>`;

        const json_data = JSON.parse(this.model.get('data'));
        const dataset = json_data as Dataset;

        const config = {
            meta: getByID(`${this.guid}-chemiscope-meta`, element),
            map: getByID(`${this.guid}-chemiscope-map`, element),
            structure: getByID(`${this.guid}-chemiscope-structure`, element),
            info: getByID(`${this.guid}-chemiscope-info`, element),
            maxStructureViewers: 4,
            // conditionally adds a settings field
            ...('settings' in json_data && { settings: json_data.settings as Partial<Settings> }),
        };

        void DefaultVisualizer.load(config, dataset)
            .then((visualizer) => {
                this.visualizer = visualizer;
            })
            .catch((e: Error) => {
                const display = getByID(`${this.guid}-error-display`, element);
                display.style.display = 'block';
                display.getElementsByTagName('p')[0].innerText = e.toString();
            });

        if (!this.model.get('has_metadata')) {
            getByID(`${this.guid}-chemiscope-meta`, element).style.display = 'none';
        }
    }

    public remove(): unknown {
        if (this.visualizer !== undefined) {
            this.visualizer.remove();
        }

        return super.remove();
    }
}

/**
 * The [[StructureView]] class renders a structure-only widget in the Jupyter
 * Notebook output window when instantiated from the Chemiscope Python package.
 */
export class StructureView extends DOMWidgetView {
    private visualizer?: StructureVisualizer;
    private guid!: string;

    public render(): void {
        this.guid = `chsp-${generateGUID()}`;

        // this function works by first rendering the widget inside `this.el`,
        // and then inserting this.el inside the HTML document.
        const element = this.el;

        addWarningHandler((message) => {
            const display = getByID(`${this.guid}-warning-display`, element);
            display.style.display = 'block';
            display.getElementsByTagName('p')[0].innerText = message;
        });

        element.innerHTML = `
        <div class="chemiscope-bootstrap">
            <div class="alert alert-warning" role="alert" id="${this.guid}-warning-display" style="display: none; font-size: 1.5em;">
                <button type="button" class="close" onclick="document.getElementById('${this.guid}-warning-display').style.display = 'none';">
                    <span aria-hidden="true">&times;</span>
                </button>
                <p></p>
            </div>
            <div class="alert alert-danger" role="alert" id="${this.guid}-error-display" style="display: none; font-size: 1.5em;">
                <button type="button" class="close" onclick="document.getElementById('${this.guid}-error-display').style.display = 'none';">
                    <span aria-hidden="true">&times;</span>
                </button>
                <p></p>
            </div>

            <div class="chemiscope-widget-one-col">
                <div class="chemiscope-structure-and-info">
                    <div class="chemiscope-meta" id="${this.guid}-chemiscope-meta"></div>
                    <div class="chemiscope-structure" id="${this.guid}-chemiscope-structure"></div>
                    <div class="chemiscope-info" id="${this.guid}-chemiscope-info"></div>
                </div>
                </div>
            </div>
        </div>`;

        const json_data = JSON.parse(this.model.get('data'));
        const dataset = json_data as Dataset;

        const config = {
            meta: getByID(`${this.guid}-chemiscope-meta`, element),
            structure: getByID(`${this.guid}-chemiscope-structure`, element),
            info: getByID(`${this.guid}-chemiscope-info`, element),
            // conditionally adds a settings field
            ...('settings' in json_data && { settings: json_data.settings as Partial<Settings> }),
        };

        // Python -> JavaScript update
        this.model.on('change:settings', this._onSettingsChanged, this);

        void StructureVisualizer.load(config, dataset)
            .then((visualizer) => {
                this.visualizer = visualizer;
                this._SaveSettings();
                this.visualizer.structure.onsettings = this._SaveSettings;
            })
            .catch((e: Error) => {
                const display = getByID(`${this.guid}-error-display`, element);
                display.style.display = 'block';
                display.getElementsByTagName('p')[0].innerText = e.toString();
            });

        if (!this.model.get('has_metadata')) {
            getByID(`${this.guid}-chemiscope-meta`, element).style.display = 'none';
        }

        // JavaScript -> Python update
        //this._emailInput.onchange = this._onInputChanged.bind(this);
    }

    private _onSettingsChanged(): void {
        console.log('onSettingsChanged');
        console.log(this.model.get('settings'));
        const settings = JSON.parse(this.model.get('settings')) as Partial<Settings>;
        console.log(settings);
        if (settings.structure !== undefined) {
            this.visualizer?.structure.applySettings(settings.structure);
        }
    }

    public _SaveSettings(): void {
        console.log('saving settings', this, this.visualizer);

        if (this.visualizer !== undefined) {
            console.log(JSON.stringify(this.visualizer.structure.saveSettings()));
            this.model.set(
                'settings',
                JSON.stringify({ structure: this.visualizer?.structure.saveSettings() })
            );
            this.model.save_changes();
        }
    }

    public remove(): unknown {
        if (this.visualizer !== undefined) {
            this.visualizer.remove();
        }

        return super.remove();
    }
}

/**
 * The [[MapView]] class renders a map-only widget in the Jupyter Notebook
 * output window when instantiated from the Chemiscope Python package.
 */
export class MapView extends DOMWidgetView {
    private visualizer?: MapVisualizer;
    private guid!: string;

    public render(): void {
        this.guid = `chsp-${generateGUID()}`;

        // this function works by first rendering the widget inside `this.el`,
        // and then inserting this.el inside the HTML document.
        const element = this.el;

        // HACK: resize Plotly when inserted. It currently render itself at full
        // width and then needs to be resized. For more on this, see
        // https://github.com/lab-cosmo/chemiscope/pull/181#discussion_r693005307
        element.addEventListener('DOMNodeInserted', () => {
            window.dispatchEvent(new Event('resize'));
        });

        addWarningHandler((message) => {
            const display = getByID(`${this.guid}-warning-display`, element);
            display.style.display = 'block';
            display.getElementsByTagName('p')[0].innerText = message;
        });

        element.innerHTML = `
        <div class="chemiscope-bootstrap">
            <div class="alert alert-warning" role="alert" id="${this.guid}-warning-display" style="display: none; font-size: 1.5em;">
                <button type="button" class="close" onclick="document.getElementById('${this.guid}-warning-display').style.display = 'none';">
                    <span aria-hidden="true">&times;</span>
                </button>
                <p></p>
            </div>
            <div class="alert alert-danger" role="alert" id="${this.guid}-error-display" style="display: none; font-size: 1.5em;">
                <button type="button" class="close" onclick="document.getElementById('${this.guid}-error-display').style.display = 'none';">
                    <span aria-hidden="true">&times;</span>
                </button>
                <p></p>
            </div>

            <div class="chemiscope-widget-one-col">
                <div class="chemiscope-meta-and-map">
                    <div class="chemiscope-meta" id="${this.guid}-chemiscope-meta"></div>
                    <div class="chemiscope-map" id="${this.guid}-chemiscope-map"></div>
                    <div class="chemiscope-info" id="${this.guid}-chemiscope-info"></div>
                </div>
                </div>
            </div>
        </div>`;

        const json_data = JSON.parse(this.model.get('data'));
        const dataset = json_data as Dataset;

        const config = {
            meta: getByID(`${this.guid}-chemiscope-meta`, element),
            map: getByID(`${this.guid}-chemiscope-map`, element),
            info: getByID(`${this.guid}-chemiscope-info`, element),
            // conditionally adds a settings field
            ...('settings' in json_data && { settings: json_data.settings as Partial<Settings> }),
        };

        void MapVisualizer.load(config, dataset)
            .then((visualizer) => {
                this.visualizer = visualizer;
            })
            .catch((e: Error) => {
                const display = getByID(`${this.guid}-error-display`, element);
                display.style.display = 'block';
                display.getElementsByTagName('p')[0].innerText = e.toString();
            });

        if (!this.model.get('has_metadata')) {
            getByID(`${this.guid}-chemiscope-meta`, element).style.display = 'none';
        }
    }

    public remove(): unknown {
        if (this.visualizer !== undefined) {
            this.visualizer.remove();
        }

        return super.remove();
    }
}
