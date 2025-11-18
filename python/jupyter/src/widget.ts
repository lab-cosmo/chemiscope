import { DOMWidgetView } from '@jupyter-widgets/base';
import Plausible from 'plausible-tracker';

import { Warnings, generateGUID, getByID } from '../../../src/utils';

// Import the CSS
import './widget.css';

import {
    DefaultConfig,
    DefaultVisualizer,
    MapVisualizer,
    StructureConfig,
    StructureVisualizer,
} from '../../../src/index';
import { Dataset, Settings, Structure, UserStructure } from '../../../src/dataset';

const PlausibleTracker = Plausible({
    domain: 'jupyter.chemiscope.org',
    // jupyter typically runs on localhost
    trackLocalhost: true,
});

interface StructureRequest {
    type: string;
    requestId: number;
    structure?: Structure | string;
    error?: string;
}

class ChemiscopeBaseView extends DOMWidgetView {
    protected visualizer?: DefaultVisualizer | StructureVisualizer | MapVisualizer;
    protected guid!: string;
    protected warnings: Warnings = new Warnings();
    protected getClassName(): string {
        return 'base-view';
    }

    // For async structure loading via Python
    private _nextRequestId = 0;
    private _pendingStructureRequests = new Map<
        number,
        { resolve: (s: Structure) => void; reject: (e: Error) => void }
    >();

    public render(): void {
        PlausibleTracker.trackPageview({
            url: (location.pathname.split('/')[1] || '') + '/' + this.getClassName(),
        });

        this.guid = `chsp-${generateGUID()}`;

        this.model.on('change:warning_timeout', () => this._updateWarningTimeout());
    }

    public remove(): unknown {
        if (this.visualizer !== undefined) {
            this.visualizer.remove();
        }

        return super.remove();
    }

    protected _bindPythonSettings(): void {
        // update settings on the JS side when they are changed in Python
        this.model.on(
            'change:settings',
            () => {
                // only trigger a visualizer update if required.
                // this is also used to avoid an infinite loop when settings are changed JS-side
                if (!this.model.get('_settings_sync')) {
                    return;
                }

                const settings = this.model.get('settings') as Partial<Settings>;

                // ignore pinned setting in jupyter, otherwise the pinned is changed
                // by JS and then overwritten the first time by Python
                delete settings.pinned;
                this.model.set('settings', settings);
                this.visualizer?.applySettings(settings);
            },
            this
        );
    }

    protected _updatePythonSettings(): void {
        if (this.visualizer !== undefined) {
            const settings = this.visualizer.saveSettings();
            // ignore pinned setting in jupyter, otherwise the pinned is changed
            // by JS and then overwritten the first time by Python
            delete settings.pinned;

            // save current settings of settings_sync
            const sync_state = this.model.get('_settings_sync') as unknown;

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

    protected _updateWarningTimeout(): void {
        const timeout = this.model.get('warning_timeout') as unknown;
        if (typeof timeout === 'number') {
            this.warnings.defaultTimeout = timeout;
        }
    }

    /**
     * Install a handler for custom messages coming from Python
     */
    protected _enableStructureMessages(): void {
        this.model.on(
            'msg:custom',
            (content: StructureRequest, _buffers: unknown[]) => {
                void _buffers; // explicitly mark as intentionally unused
                if (!content || typeof content !== 'object') {
                    return;
                }

                if (content.type === 'load-structure-result') {
                    const requestId = content.requestId;
                    const pending = this._pendingStructureRequests.get(requestId);
                    if (pending) {
                        this._pendingStructureRequests.delete(requestId);
                        // python can send either a string (JSON-serialized) or an object
                        const raw = content.structure;
                        const structure: Structure =
                            typeof raw === 'string'
                                ? (JSON.parse(raw) as Structure)
                                : (raw as Structure);
                        pending.resolve(structure);
                    }
                } else if (content.type === 'load-structure-error') {
                    const requestId = content.requestId;
                    const pending = this._pendingStructureRequests.get(requestId);
                    if (pending) {
                        this._pendingStructureRequests.delete(requestId);
                        pending.reject(
                            new Error(content.error || 'unknown error while loading structure')
                        );
                    }
                    if (content.error) {
                        this.warnings.sendMessage(`Error loading structure: ${content.error}`);
                    }
                }
            },
            this
        );
    }

    /**
     * Ask Python to load one structure given the user-defined `data`.
     * `data` is typically the filename or an object containing it.
     */
    protected _requestStructureFromPython(
        index: number,
        structure: UserStructure
    ): Promise<Structure> {
        const requestId = this._nextRequestId++;
        return new Promise<Structure>((resolve, reject) => {
            this._pendingStructureRequests.set(requestId, { resolve, reject });

            this.model.send({
                type: 'load-structure',
                requestId,
                index,
                data: structure.data as string,
            });
        });
    }

    /**
     * If the dataset uses `UserStructure`, wrap the structure config so that
     * chemiscope will ask Python to load each structure on demand.
     */
    protected _attachStructureLoaderToConfig(
        config: DefaultConfig | StructureConfig,
        data: Dataset
    ): void {
        const structures = data.structures;
        if (!Array.isArray(structures) || structures.length === 0) {
            return;
        }

        const structure = structures[0];
        const hasUserData = structure && typeof structure === 'object' && 'data' in structure;
        if (!hasUserData) {
            // Normal (already-expanded) structures, nothing to do
            return;
        }

        // We are in the dynamic UserStructure case.
        // Make sure we listen for Python replies
        this._enableStructureMessages();

        // StructureConfig.loadStructure gets the full UserStructure as second arg
        config.loadStructure = async (index: number, raw: unknown): Promise<Structure> => {
            const wrapper = raw as UserStructure;

            const promise = this._requestStructureFromPython(index, wrapper);
            return promise;
        };
    }
}

/**
 * The {@link ChemiscopeView} class renders the Chemiscope App as a widget in the
 * Jupyter Notebook output window when instantiated from the Chemiscope Python
 * package.
 */
export class ChemiscopeView extends ChemiscopeBaseView {
    protected visualizer?: DefaultVisualizer;
    protected getClassName(): string {
        return 'chemiscope-view';
    }

    public render(): void {
        super.render();

        // this function works by first rendering the widget inside `this.el`,
        // and then inserting this.el inside the HTML document.
        const element = this.el;

        this._updateWarningTimeout();
        this.warnings.addHandler((message, timeout?) => {
            displayWarning(message, element, this.guid, timeout);
        });

        element.innerHTML = `
        <div>
            <div class="alert alert-warning alert-dismissible pop-on-top" role="alert" id="${this.guid}-warning-display" style="display: none; font-size: 1.5em;">
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

            <div class="chemiscope-viewer-two-col">
                <div class="chemiscope-meta-and-map">
                    <div class="chemiscope-meta" id="${this.guid}-chemiscope-meta"></div>
                    <div class="chemiscope-map" id="${this.guid}-chemiscope-map"></div>
                </div>
                <div class="chemiscope-structure-and-info">
                    <div class="chemiscope-structure" id="${this.guid}-chemiscope-structure"></div>
                    <div class="chemiscope-info" id="${this.guid}-chemiscope-info"></div>
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

        const data = parseJsonWithNaN(this.model.get('value') as string) as Dataset;

        // If `data.structures` is an array of `UserStructure`, this
        // will wrap the structure config with a loader that calls Python.
        this._attachStructureLoaderToConfig(config, data);

        void DefaultVisualizer.load(config, data, this.warnings)
            .then((visualizer) => {
                this.visualizer = visualizer;
                // update the Python side settings whenever a setting changes
                this.visualizer.onSettingChange(() => this._updatePythonSettings());
                // and set them to the initial value right now
                this._updatePythonSettings();
            })
            // eslint-disable-next-line @typescript-eslint/use-unknown-in-catch-callback-variable
            .catch((e: Error) => {
                // eslint-disable-next-line no-console
                console.error(e);

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
 * The {@link StructureView} class renders a structure-only widget in the Jupyter
 * Notebook output window when instantiated from the Chemiscope Python package.
 */
export class StructureView extends ChemiscopeBaseView {
    protected visualizer?: StructureVisualizer;
    protected getClassName(): string {
        return 'structure-view';
    }

    public render(): void {
        super.render();

        // this function works by first rendering the widget inside `this.el`,
        // and then inserting this.el inside the HTML document.
        const element = this.el;

        this._updateWarningTimeout();
        this.warnings.addHandler((message, timeout?) => {
            displayWarning(message, element, this.guid, timeout);
        });

        element.innerHTML = `
        <div>
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

            <div class="chemiscope-viewer-one-col">
                <div class="chemiscope-structure-and-info">
                    <div class="chemiscope-meta" id="${this.guid}-chemiscope-meta"></div>
                    <div class="chemiscope-structure" id="${this.guid}-chemiscope-structure"></div>
                    <div class="chemiscope-info" id="${this.guid}-chemiscope-info"></div>
                </div>
            </div>
        </div>`;

        const config = {
            meta: getByID(`${this.guid}-chemiscope-meta`, element),
            structure: getByID(`${this.guid}-chemiscope-structure`, element),
            info: getByID(`${this.guid}-chemiscope-info`, element),
        };

        this._bindPythonSettings();

        const data = parseJsonWithNaN(this.model.get('value') as string) as Dataset;

        // If `data.structures` is an array of `UserStructure`, this
        // will wrap the structure config with a loader that calls Python.
        this._attachStructureLoaderToConfig(config, data);

        void StructureVisualizer.load(config, data, this.warnings)
            .then((visualizer) => {
                this.visualizer = visualizer;

                // update the Python side settings whenever a setting changes
                this.visualizer.onSettingChange(() => this._updatePythonSettings());
                // and set them to the initial value right now
                this._updatePythonSettings();
            })
            // eslint-disable-next-line @typescript-eslint/use-unknown-in-catch-callback-variable
            .catch((e: Error) => {
                // eslint-disable-next-line no-console
                console.error(e);

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
 * The {@link MapView} class renders a map-only widget in the Jupyter Notebook
 * output window when instantiated from the Chemiscope Python package.
 */
export class MapView extends ChemiscopeBaseView {
    protected visualizer?: MapVisualizer;
    protected getClassName(): string {
        return 'map-view';
    }

    public render(): void {
        super.render();

        // this function works by first rendering the widget inside `this.el`,
        // and then inserting this.el inside the HTML document.
        const element = this.el;

        this._updateWarningTimeout();
        this.warnings.addHandler((message, timeout?) => {
            displayWarning(message, element, this.guid, timeout);
        });

        element.innerHTML = `
        <div>
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

            <div class="chemiscope-viewer-one-col">
                <div class="chemiscope-meta-and-map">
                    <div class="chemiscope-meta" id="${this.guid}-chemiscope-meta"></div>
                    <div class="chemiscope-map" id="${this.guid}-chemiscope-map"></div>
                    <div class="chemiscope-info" id="${this.guid}-chemiscope-info"></div>
                </div>
            </div>
        </div>`;

        const config = {
            meta: getByID(`${this.guid}-chemiscope-meta`, element),
            map: getByID(`${this.guid}-chemiscope-map`, element),
            info: getByID(`${this.guid}-chemiscope-info`, element),
        };

        this._bindPythonSettings();

        const data = parseJsonWithNaN(this.model.get('value') as string) as Dataset;
        void MapVisualizer.load(config, data, this.warnings)
            .then((visualizer) => {
                this.visualizer = visualizer;

                // update the Python side settings whenever a setting changes
                this.visualizer.onSettingChange(() => this._updatePythonSettings());
                // and set them to the initial value right now
                this._updatePythonSettings();
            })
            // eslint-disable-next-line @typescript-eslint/use-unknown-in-catch-callback-variable
            .catch((e: Error) => {
                // eslint-disable-next-line no-console
                console.error(e);

                const display = getByID(`${this.guid}-error-display`, element);
                display.style.display = 'block';
                display.getElementsByTagName('p')[0].innerText = e.toString();
            });

        if (!this.model.get('has_metadata')) {
            getByID(`${this.guid}-chemiscope-meta`, element).style.display = 'none';
        }
    }
}

function displayWarning(
    message: string,
    element: HTMLElement,
    guid: string,
    timeout: number = 4000
) {
    if (timeout < 0) return;

    const display = getByID(`${guid}-warning-display`, element);
    display.style.display = 'block';
    display.getElementsByTagName('p')[0].innerText = message;

    if (timeout > 0) {
        // automatically remove the warning after a set timeout
        setTimeout(() => {
            display.style.display = 'none';
        }, timeout);
    }
}

/**
 * Allow NaN in the JSON file. They are not part of the spec, but Python's json
 * module output them, and they can be useful.
 */
function parseJsonWithNaN(text: string): unknown {
    return JSON.parse(
        text.replace(/\bNaN\b/g, '"***NaN***"'),
        (key: string, value: unknown): unknown => {
            return value === '***NaN***' ? NaN : value;
        }
    ) as unknown;
}
