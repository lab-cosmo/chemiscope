import { DOMWidgetView } from '@jupyter-widgets/base';
import { addWarningHandler, generateGUID, getByID } from '../../../src/utils';

// Import the CSS
import './widget.css';
import './chemiscope-bootstrap.less';
import 'bootstrap/dist/js/bootstrap.min.js';

import { DefaultVisualizer, MapVisualizer, StructureVisualizer } from '../../../src/index';
import { Dataset, Settings } from '../../../src/dataset';


class ChemiscopeBaseView extends DOMWidgetView {
    protected visualizer?: any;
    protected guid!: string;
    
    public remove(): unknown {
        if (this.visualizer !== undefined) {
            this.visualizer.remove();
        }

        return super.remove();
    }

    protected _bindPythonSettings(): void {
        // update settings on the JS side when they are changed in Python
        this.model.on('change:settings', () => {
            // only trigger a visualizer update if required.
            // this is also used to avoid an infinite loop when settings are changed JS-side
            if (!this.model.get('_settings_sync')) { 
                return;
            }

            const settings = this.model.get('settings') as Partial<Settings>;
            
            // ignore pinned setting in jupyter, otherwise the pinned is changed
            // by JS and then overwritten the first time by Python
            delete settings.pinned;
            this.model.set('settings', settings)
            this.visualizer?.applySettings(settings);
        }, this);
    }

    protected _updatePythonSettings(): void {
        if (this.visualizer !== undefined) {            
            const settings = this.visualizer.saveSettings();
            // ignore pinned setting in jupyter, otherwise the pinned is changed
            // by JS and then overwritten the first time by Python
            delete settings.pinned;
            
            // save current settings of settings_sync
            let sync_state = this.model.get('_settings_sync'); 

            // signals that updating the Python state shouldn't trigger a re-update.
            // this is a workaround because it seems that settings:change doesn't know
            // if it's triggered from JS or from Python, so we need an extra flag to avoid a loop
            this.model.set('_settings_sync', false);  
            this.model.save_changes(); 
            this.model.set('settings', settings);
            this.model.save_changes();
            this.model.set('_settings_sync', sync_state);       
            this.model.save_changes();
        }
    }
}

/**
 * The [[ChemiscopeView]] class renders the Chemiscope App as a widget in the
 * Jupyter Notebook output window when instantiated from the Chemiscope Python
 * package.
 */
export class ChemiscopeView extends CSBaseView {
    protected visualizer?: DefaultVisualizer;
    
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

        const config = {
            meta: getByID(`${this.guid}-chemiscope-meta`, element),
            map: getByID(`${this.guid}-chemiscope-map`, element),
            structure: getByID(`${this.guid}-chemiscope-structure`, element),
            info: getByID(`${this.guid}-chemiscope-info`, element),
            maxStructureViewers: 4,
        };

        this._bindPythonSettings();

        const data = JSON.parse(this.model.get('data') as string) as Dataset;
        void DefaultVisualizer.load(config, data)
            .then((visualizer) => {
                this.visualizer = visualizer;
                // update the Python side settings whenever a setting changes
                this.visualizer.onSettingChange(() => this._updatePythonSettings());
                // and set them to the initial value right now
                this._updatePythonSettings();
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
}

/**
 * The [[StructureView]] class renders a structure-only widget in the Jupyter
 * Notebook output window when instantiated from the Chemiscope Python package.
 */
export class StructureView extends CSBaseView {
    protected visualizer?: StructureVisualizer;
    
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

        const config = {
            meta: getByID(`${this.guid}-chemiscope-meta`, element),
            structure: getByID(`${this.guid}-chemiscope-structure`, element),
            info: getByID(`${this.guid}-chemiscope-info`, element),
        };

        this._bindPythonSettings();

        const data = JSON.parse(this.model.get('data') as string) as Dataset;        
        void StructureVisualizer.load(config, data)
            .then((visualizer) => {
                this.visualizer = visualizer;

                // update the Python side settings whenever a setting changes
                this.visualizer.onSettingChange(() => this._updatePythonSettings());
                // and set them to the initial value right now
                this._updatePythonSettings();
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
}

/**
 * The [[MapView]] class renders a map-only widget in the Jupyter Notebook
 * output window when instantiated from the Chemiscope Python package.
 */
export class MapView extends CSBaseView {
    protected visualizer?: MapVisualizer;
    
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

        const config = {
            meta: getByID(`${this.guid}-chemiscope-meta`, element),
            map: getByID(`${this.guid}-chemiscope-map`, element),
            info: getByID(`${this.guid}-chemiscope-info`, element),
        };

        this._bindPythonSettings();

        const data = JSON.parse(this.model.get('data') as string) as Dataset;
        void MapVisualizer.load(config, data)
            .then((visualizer) => {
                this.visualizer = visualizer;

                // update the Python side settings whenever a setting changes
                this.visualizer.onSettingChange(() => this._updatePythonSettings());
                // and set them to the initial value right now
                this._updatePythonSettings();
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
}
