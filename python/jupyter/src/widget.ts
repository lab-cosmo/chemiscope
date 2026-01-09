import { DOMWidgetView } from '@jupyter-widgets/base';
import { JSONValue } from '@lumino/coreutils';
import Plausible from 'plausible-tracker';

import { Warnings, binarySearch, generateGUID, getByID } from '../../../src/utils';

// Import the CSS
import './widget.css';

import {
    DefaultConfig,
    DefaultVisualizer,
    DisplayTarget,
    Indexes,
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
    data?: string;
    error?: string;
}

interface ScreenshotRequest {
    type: string;
    requestId: number;
    target: 'map' | 'structure';
}

interface StructureSequenceRequest {
    type: string;
    requestId: number;
    indices: number[];
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

    // Flag to prevent infinite loops when updating settings from Python
    protected _updatingFromPython = false;

    public render(): void {
        PlausibleTracker.trackPageview({
            url: (location.pathname.split('/')[1] || '') + '/' + this.getClassName(),
        });

        this.guid = `chsp-${generateGUID()}`;

        this.model.on('change:warning_timeout', () => this._updateWarningTimeout());
        this._enableScreenshotMessages();
    }

    public remove(): unknown {
        if (this.visualizer !== undefined) {
            this.visualizer.remove();
        }

        return super.remove();
    }

    protected _initializeVisualizer(
        visualizer: DefaultVisualizer | StructureVisualizer | MapVisualizer
    ): void {
        this.visualizer = visualizer;

        const settings = this.model.get('settings') as Partial<Settings>;
        if (settings) {
            this.visualizer.applySettings(settings);
        }

        // update the Python side settings whenever a setting changes
        this.visualizer.onSettingChange(() => {
            if (!this._updatingFromPython) {
                this._updatePythonSettings();
            }
        });
        // and set them to the initial value right now
        this._updatePythonSettings();
        this._bindSelection();
    }

    protected _bindPythonSettings(): void {
        // update settings on the JS side when they are changed in Python
        this.model.on(
            'change:settings',
            () => {
                // only trigger a visualizer update if required.
                // this is also used to avoid an infinite loop when settings are changed JS-side
                if (this._updatingFromPython) {
                    return;
                }

                const settingsRef = this.model.get('settings') as Partial<Settings>;
                const settings = { ...settingsRef };

                // Handle pinned: only apply if different from current state.
                // This prevents the destructive reset loop while allowing explicit updates.
                if (
                    this.visualizer &&
                    'structure' in this.visualizer &&
                    Array.isArray(settings.pinned)
                ) {
                    const target = this.visualizer.saveSettings().target as DisplayTarget;
                    const currentPinned = this.visualizer.structure
                        .pinned()
                        .map((value) => (target === 'atom' ? value.environment : value.structure));

                    const pinned = settings.pinned as number[];
                    const pinnedChanged =
                        pinned.length !== currentPinned.length ||
                        pinned.some((val, idx) => val !== currentPinned[idx]);

                    if (!pinnedChanged) {
                        delete settings.pinned;
                    }
                } else {
                    delete settings.pinned;
                }

                this._updatingFromPython = true;
                try {
                    this.visualizer?.applySettings(settings);
                    // sync back the full settings to Python, so that they are available
                    // for saving
                    this._updatePythonSettings();

                    // also sync back the selection, which might have changed if pinned changed
                    if (this.visualizer) {
                        this._updatePythonSelection(this.visualizer.info.indexes, true);
                    }
                } catch (e) {
                    this.warnings.sendMessage(`Error setting state: ${e}`);
                } finally {
                    this._updatingFromPython = false;
                }
            },
            this
        );
    }

    protected _updatePythonSelection(indexes: Indexes, force = false): void {
        if (this._updatingFromPython && !force) {
            return;
        }

        const wasUpdating = this._updatingFromPython;
        this._updatingFromPython = true;
        try {
            const currentSelected = this.model.get('selected_ids') as
                | {
                      structure: number;
                      atom?: number;
                  }
                | undefined;

            const selectedChanged =
                !currentSelected ||
                currentSelected.structure !== indexes.structure ||
                currentSelected.atom !== indexes.atom;

            if (selectedChanged) {
                this.model.set('selected_ids', {
                    structure: indexes.structure,
                    atom: indexes.atom,
                });
            }

            if (this.visualizer && 'structure' in this.visualizer) {
                const activeViewer = this.visualizer.structure.activeIndex;
                if (this.model.get('active_viewer') !== activeViewer) {
                    this.model.set('active_viewer', activeViewer);
                }
            }

            this._updatePythonSettings();
            this.model.save_changes();
        } finally {
            this._updatingFromPython = wasUpdating;
        }
    }

    protected _bindSelection(): void {
        if (!this.visualizer) {
            return;
        }

        // JS -> Python
        const updatePython = (indexes: Indexes) => this._updatePythonSelection(indexes);

        if ('structure' in this.visualizer) {
            const originalOnSelect = this.visualizer.structure.onselect;
            this.visualizer.structure.onselect = (indexes) => {
                if (originalOnSelect) {
                    originalOnSelect(indexes);
                }
                updatePython(indexes);
            };

            const originalActiveChanged = this.visualizer.structure.activeChanged;
            this.visualizer.structure.activeChanged = (guid, indexes) => {
                if (originalActiveChanged) {
                    originalActiveChanged(guid, indexes);
                }
                updatePython(indexes);
            };

            const originalRemoveViewer = this.visualizer.structure.removeViewer.bind(
                this.visualizer.structure
            );
            this.visualizer.structure.removeViewer = (guid) => {
                originalRemoveViewer(guid);
                if (this.visualizer) {
                    updatePython(this.visualizer.info.indexes);
                }
            };
        }

        if ('map' in this.visualizer) {
            const originalOnSelect = this.visualizer.map.onselect;
            this.visualizer.map.onselect = (indexes) => {
                if (originalOnSelect) {
                    originalOnSelect(indexes);
                }
                updatePython(indexes);
            };

            const originalActiveChanged = this.visualizer.map.activeChanged;
            this.visualizer.map.activeChanged = (guid, indexes) => {
                if (originalActiveChanged) {
                    originalActiveChanged(guid, indexes);
                }
                updatePython(indexes);
            };
        }

        if (this.visualizer.info) {
            const originalOnChange = this.visualizer.info.onchange;
            this.visualizer.info.onchange = (indexes) => {
                if (originalOnChange) {
                    originalOnChange(indexes);
                }
                updatePython(indexes);
            };
        }

        // Python -> JS
        this.model.on(
            'change:selected_ids',
            () => {
                if (!this.visualizer) {
                    return;
                }

                if (this._updatingFromPython) {
                    return;
                }

                const selected = this.model.get('selected_ids') as
                    | {
                          structure?: number;
                          atom?: number;
                      }
                    | undefined;

                if (!selected) {
                    return;
                }

                this._updatingFromPython = true;
                try {
                    const currentIndexes = this.visualizer.info.indexes;
                    let structure = selected.structure;
                    let atom = selected.atom;

                    if (structure === undefined) {
                        structure = currentIndexes.structure;
                    }

                    if (atom === undefined) {
                        atom = currentIndexes.atom;
                    }

                    // Validate atom index

                    if (atom !== undefined) {
                        const activeAtoms = this.visualizer.indexer.activeAtoms(structure);

                        if (activeAtoms && activeAtoms.length > 0) {
                            if (binarySearch(activeAtoms, atom) === -1) {
                                // Reset to 0 if valid, else first active
                                if (binarySearch(activeAtoms, 0) !== -1) {
                                    atom = 0;
                                } else {
                                    atom = activeAtoms[0];
                                }
                            }
                        }
                    }

                    const target = this.visualizer.saveSettings().target as DisplayTarget;
                    const indexes = this.visualizer.indexer.fromStructureAtom(
                        target,
                        structure,
                        atom
                    );

                    if (indexes !== undefined) {
                        this.visualizer.select(indexes);
                    } else {
                        this.warnings.sendMessage(
                            `Invalid selection request: structure=${selected.structure}, atom=${selected.atom}`
                        );
                    }
                } finally {
                    this._updatingFromPython = false;
                }
            },
            this
        );

        // Python -> JS (active viewer)
        this.model.on(
            'change:active_viewer',
            () => {
                if (!this.visualizer || !('structure' in this.visualizer)) {
                    return;
                }

                if (this._updatingFromPython) {
                    return;
                }

                const activeViewer = this.model.get('active_viewer') as number;

                this._updatingFromPython = true;
                try {
                    this.visualizer.structure.activeIndex = activeViewer;
                    if ('map' in this.visualizer) {
                        this.visualizer.map.setActive(this.visualizer.structure.active);
                    }
                } finally {
                    this._updatingFromPython = false;
                }
            },
            this
        );

        // Initialize Python state if needed
        const currentPython = this.model.get('selected_ids') as {
            structure?: number;
        };

        if (!currentPython || currentPython.structure === undefined) {
            // we use the info panel as the source of truth for the current selection
            const indexes = this.visualizer.info.indexes;
            updatePython(indexes);
        }
    }

    protected _updatePythonSettings(): void {
        if (this.visualizer !== undefined) {
            const settings = this.visualizer.saveSettings();
            // ignore pinned setting in jupyter, otherwise the pinned is changed
            // by JS and then overwritten the first time by Python
            delete settings.pinned;

            this.model.set('settings', settings);
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
    protected _enableScreenshotMessages(): void {
        this.model.on(
            'msg:custom',
            (content: ScreenshotRequest | StructureSequenceRequest, _buffers: unknown[]) => {
                void _buffers;
                if (!content) {
                    return;
                }

                if (content.type === 'save-image') {
                    const req = content as ScreenshotRequest;
                    this._handleSaveImage(req);
                } else if (content.type === 'save-structure-sequence') {
                    const req = content as StructureSequenceRequest;
                    void this._handleStructureSequence(req);
                }
            },
            this
        );
    }

    private _handleSaveImage(content: ScreenshotRequest) {
        const requestId = content.requestId;
        const target = content.target;

        const sendResult = (data: string) => {
            this.model.send({
                type: 'save-image-result',
                requestId: requestId,
                data: data,
            });
        };

        const sendError = (error: string) => {
            this.model.send({
                type: 'save-image-error',
                requestId: requestId,
                error: error,
            });
        };

        if (!this.visualizer) {
            sendError('visualizer not loaded');
            return;
        }

        if (target === 'map') {
            if ('map' in this.visualizer) {
                this.visualizer.map
                    .exportPNG()
                    .then(sendResult)
                    .catch((e: unknown) => sendError(`${e}`));
            } else {
                sendError('no map in this visualizer');
            }
        } else if (target === 'structure') {
            if ('structure' in this.visualizer) {
                try {
                    const data = this.visualizer.structure.exportActivePNG();
                    sendResult(data);
                } catch (e) {
                    sendError((e as Error).toString());
                }
            } else {
                sendError('no structure in this visualizer');
            }
        } else {
            sendError(`unknown target ${target}`);
        }
    }

    private async _handleStructureSequence(content: StructureSequenceRequest) {
        const requestId = content.requestId;
        const indices = content.indices;

        if (!this.visualizer || !('structure' in this.visualizer)) {
            this.model.send({
                type: 'save-structure-sequence-error',
                requestId: requestId,
                error: 'no structure visualizer available',
            });
            return;
        }

        const indexer = this.visualizer.indexer;
        const structureViewer = this.visualizer.structure;

        for (const index of indices) {
            try {
                const target = this.visualizer.saveSettings().target as DisplayTarget;
                const indexes = indexer.fromStructure(index, target);

                if (indexes) {
                    await structureViewer.show(indexes);
                    await new Promise((r) => requestAnimationFrame(r));

                    const data = structureViewer.exportActivePNG();
                    this.model.send({
                        type: 'save-structure-sequence-result',
                        requestId: requestId,
                        index: index,
                        data: data,
                    });
                }
            } catch (e) {
                // eslint-disable-next-line no-console
                console.error(`Failed to save frame ${index}`, e);
            }
        }

        this.model.send({
            type: 'save-structure-sequence-done',
            requestId: requestId,
        });
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
            let limit = 16;

            // Limit concurrency during playback to prevent flooding

            if (this.visualizer && 'info' in this.visualizer && this.visualizer.info.isPlaying) {
                limit = 1;
            }

            // avoid piling up too many requests during traj. playback

            if (this._pendingStructureRequests.size >= limit) {
                // eslint-disable-next-line no-console
                console.warn(
                    `Skipping structure ${index} - ${structure.data}. Increase playback delay.`
                );
            } else {
                // queue a request for the structure
                this._pendingStructureRequests.set(requestId, { resolve, reject });

                this.model.send({
                    type: 'load-structure',
                    requestId,
                    index,
                    data: structure.data as JSONValue,
                });
            }
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
                this._initializeVisualizer(visualizer);
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
                this._initializeVisualizer(visualizer);
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
                this._initializeVisualizer(visualizer);
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
