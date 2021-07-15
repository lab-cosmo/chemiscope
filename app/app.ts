// load bootstrap: this needs to come first so that CSS files are loaded in the
// right order
require('bootstrap');
require('bootstrap/dist/css/bootstrap.min.css');

import { getByID, addWarningHandler } from '../src/utils';
import { Dataset, Structure } from '../src/dataset';
import { version, DefaultVisualizer, Settings } from '../src/index';

import { inflate } from 'pako';

// load CSS for the app
require('./app.css');

interface Configuration {
    /// saved settings to apply to the visualizer
    settings: Partial<Settings>;
    /// optional callback to load the structures on demand.
    loadStructure?: (index: number, structure: unknown) => Structure;
}

function cleanOutdatedMessages() {
    const displayWarning = getByID('warning-display');
    displayWarning.style.display = 'none';
    displayWarning.getElementsByTagName('p')[0].innerText = '';

    const displayError = document.getElementById('error-display');
    if (displayError != null) {
        displayError.style.display = 'none';
        displayError.getElementsByTagName('p')[0].innerText = '';
        const stacktrace = displayError.getElementsByTagName('details')[0];
        stacktrace.getElementsByTagName('p')[0].innerText = '';
    }
}

export class ChemiscopeApp {
    /// Instance of the chemiscope visualizer
    private visualizer?: DefaultVisualizer;
    /// Path/URL of the dataset currently displayed
    private dataset?: string;
    /// CSS style sheet to hide the setting in loading panel about on-demand
    /// loading
    private hideOnDemandStructures: HTMLStyleElement;

    /**
     * Create a new instance of the chemiscope application.
     */
    constructor(id: string) {
        // show the version of chemiscope currently running
        const versionDisplay = getByID('chemiscope-version');
        versionDisplay.innerText = `version ${version()}`;

        // handle warnings
        addWarningHandler((message) => displayWarning(message));

        // when the window is resized, change the size available to the info
        // widget
        window.addEventListener('resize', updateInfoWidgetHeight);

        // setup the main HTML
        const root = getByID(id);
        root.innerHTML = `<main class="container-fluid">
            <div class="row">
                <div class="col-md-7" style="padding: 0">
                    <div class="embed-responsive embed-responsive-1by1">
                        <div id="chemiscope-meta"></div>
                        <div id="chemiscope-map" class="embed-responsive-item" style="position: absolute"></div>
                    </div>
                </div>

                <div class="col-md-5" style="padding: 0">
                    <div class="embed-responsive embed-responsive-5by7">
                        <div class="embed-responsive-item">
                            <!-- height: 0 below is a hack to force safari to
                            respect height: 100% on the children
                            https://github.com/philipwalton/flexbugs/issues/197#issuecomment-378908438
                            -->
                            <div id="chemiscope-structure" style="height: 0"></div>
                            <div id="chemiscope-info"></div>
                        </div>
                    </div>
                </div>
            </div>
        </main>`;

        this._setupLoadSaveMenu();

        this.hideOnDemandStructures = document.createElement('style');
        document.head.appendChild(this.hideOnDemandStructures);
        this.hideOnDemandStructures.sheet!.insertRule('.hide-on-demand-structures {display: none}');
        this.hideOnDemandStructures.sheet!.disabled = false;
    }

    /**
     * Load the example input with the given name, and display it
     */
    public async loadExample(example: string): Promise<void> {
        startLoading();
        this.dataset = `/examples/${example}.json.gz`;

        const response = await fetch(this.dataset);
        if (!response.ok) {
            throw Error(`unable to load file at ${this.dataset}`);
        }
        const buffer = await response.arrayBuffer();
        const dataset = readJSON(this.dataset, buffer);

        let config: Configuration = {
            settings: dataset.settings || {},
            loadStructure: undefined,
        };
        if (example == 'Azaphenacenes') {
            // example of asynchronous structure loading
            config.loadStructure = (_: number, structure: any) => {
                return JSON.parse(
                    $.ajax({
                        type: 'GET',
                        url: `examples/${structure.data}`,
                        // this is getting deprecated, but the best option for now
                        async: false,
                    }).responseText
                );
            };
        }

        await this.load(config, dataset);
    }

    /**
     * Load the given `dataset` with the specified `configuration`.
     */

    public async load(configuration: Configuration, dataset: Dataset): Promise<void> {
        cleanOutdatedMessages();
        const config = {
            map: 'chemiscope-map',
            info: 'chemiscope-info',
            meta: 'chemiscope-meta',
            structure: 'chemiscope-structure',
            settings: configuration.settings,
            loadStructure: configuration.loadStructure,
        };

        // show/hide setting related to on-demand structure loading
        this.hideOnDemandStructures.sheet!.disabled = configuration.loadStructure !== undefined;

        if (this.visualizer !== undefined) {
            this.visualizer.remove();
        }

        this.visualizer = await DefaultVisualizer.load(config, dataset);

        this.visualizer.structure.positionSettingsModal = (rect) => {
            const structureRect = getByID('chemiscope-structure').getBoundingClientRect();
            return {
                top: structureRect.top,
                left: structureRect.left - rect.width - 25,
            };
        };

        this.visualizer.map.positionSettingsModal = (rect) => {
            const mapRect = getByID('chemiscope-map').getBoundingClientRect();

            let left;
            if (window.innerWidth < 1400) {
                // clip modal to the right if it overflows
                left = window.innerWidth - rect.width - 10;
            } else {
                left = mapRect.left + mapRect.width + 25;
            }

            return {
                top: mapRect.top,
                left: left,
            };
        };

        updateInfoWidgetHeight();
        stopLoading();
    }

    /**
     * Setup all callbacks & style related to the Load/Save menu
     */
    private _setupLoadSaveMenu(): void {
        // Loading new dataset
        const loadDataset = getByID<HTMLInputElement>('load-dataset');
        const loadSaveModal = getByID('load-save');
        const closeLoadSaveModal = getByID('close-load-save-modal');
        loadDataset.onchange = () => {
            // remove closing animation on the modal to close it with JS
            loadSaveModal.classList.remove('fade');
            closeLoadSaveModal.click();
            startLoading();
            const file = loadDataset.files![0];
            this.dataset = file.name;
            readFile(file, (result) => {
                const dataset = readJSON(file.name, result);
                const config = {
                    settings: dataset.settings || {},
                };
                this.load(config, dataset);
                // clear the selected file name to make sure 'onchange' is
                // called again if the user loads a file a the same path
                // multiple time
                loadDataset.value = '';
                loadSaveModal.classList.add('fade');
            });
        };
        // Saving the current dataset
        const saveDataset = getByID('save-dataset');
        const saveDatasetName = getByID<HTMLInputElement>('save-dataset-name');
        const includeSettings = getByID<HTMLInputElement>('save-dataset-settings');
        const includeStructures = getByID<HTMLInputElement>('save-dataset-structures');
        saveDataset.onclick = () => {
            if (this.visualizer === undefined) {
                return;
            }

            const dataset: any = this.visualizer.dataset(includeStructures.checked);
            if (includeSettings.checked) {
                dataset.settings = this.visualizer.saveSettings();
            }
            startDownload(saveDatasetName.value, stringifyJsonWithNaN(dataset));
            closeLoadSaveModal.click();
        };
        // loading saved settings
        const loadSettings = getByID<HTMLInputElement>('load-settings');
        loadSettings.onchange = () => {
            loadSaveModal.classList.remove('fade');
            closeLoadSaveModal.click();
            startLoading();
            const file = loadSettings.files![0];
            readFile(file, (result) => {
                if (this.visualizer === undefined) {
                    return;
                }

                this.visualizer.applySettings(readJSON(file.name, result));
                // clear the selected file name to make sure 'onchange' is
                // called again if the user loads a file a the same path
                // multiple time
                loadSettings.value = '';
                stopLoading();
                loadSaveModal.classList.add('fade');
            });
        };

        // Saving the current settings values
        const saveSettings = getByID('save-settings');
        const saveSettingsName = getByID<HTMLInputElement>('save-settings-name');
        const saveSettingsMap = getByID<HTMLInputElement>('save-settings-map');
        const saveSettingsStructure = getByID<HTMLInputElement>('save-settings-structure');
        const saveSettingsSelected = getByID<HTMLInputElement>('save-settings-selected');
        saveSettings.onclick = () => {
            if (this.visualizer === undefined) {
                return;
            }

            const settings: Partial<Settings> = this.visualizer.saveSettings();
            if (!saveSettingsMap.checked) {
                delete settings.map;
            }
            if (!saveSettingsStructure.checked) {
                delete settings.structure;
            }
            if (!saveSettingsSelected.checked) {
                delete settings.pinned;
            }
            startDownload(saveSettingsName.value, JSON.stringify(settings));
            closeLoadSaveModal.click();
        };
    }
}

function displayWarning(message: string) {
    const display = getByID('warning-display');
    display.style.display = 'block';
    display.getElementsByTagName('p')[0].innerText = message;
}

function startLoading() {
    getByID('loading').style.display = 'block';

    const main = document.getElementsByTagName('main')[0];
    main.onclick = () => {};
    main.style.opacity = '0.3';
}

function stopLoading() {
    getByID('loading').style.display = 'none';

    const main = document.getElementsByTagName('main')[0];
    main.onclick = null;
    main.style.opacity = '1';
}

/** Read JSON or gzipped JSON and return the parsed object */
function readJSON(path: string, buffer: ArrayBuffer): any {
    let text;
    if (path.endsWith('.gz')) {
        text = inflate(new Uint8Array(buffer), { to: 'string' });
    } else {
        const decoder = new TextDecoder('utf-8');
        text = decoder.decode(buffer);
    }
    return parseJsonWithNaN(text);
}

/**
 * Allow NaN in the JSON file. They are not part of the spec, but Python's json
 * module output them, and they can be useful.
 */
function parseJsonWithNaN(text: string): any {
    return JSON.parse(text.replace(/\bNaN\b/g, '"***NaN***"'), (key, value) => {
        return value === '***NaN***' ? NaN : value;
    });
}

/** Write NaN in the JSON file. */
function stringifyJsonWithNaN(object: any): string {
    const string = JSON.stringify(object, (key, value) => {
        return typeof value === 'number' && isNaN(value) ? '***NaN***' : value;
    });
    return string.replace(/"\*\*\*NaN\*\*\*"/g, 'NaN');
}

/**
 * Read a `file` with FileReader, and use `callback` once loaded.
 *
 * The callback should take a single ArrayBuffer parameter.
 *
 * @param  {File}     file     user-provided file to read
 * @param  {Function} callback callback to use once the file is read
 */
function readFile(file: File, callback: (data: ArrayBuffer) => void): void {
    const reader = new FileReader();
    reader.onload = () => {
        if (reader.error) {
            throw Error(`could not read ${file.name}: ${reader.error}`);
        }
        if (reader.result) {
            callback(reader.result as ArrayBuffer);
        }
    };
    reader.readAsArrayBuffer(file);
}

/**
 * Ensure that the height of the property table is set up in a compatible way
 * with the use of embed-responsive. embed-responsive is used to ensure no
 * vertical or horizontal scroll bar appear on the whole page while keeping
 * know aspect ratio for the main panels.
 */
function updateInfoWidgetHeight(): void {
    const height = getByID('chemiscope-structure').getBoundingClientRect().height;

    const infoTables = document.getElementsByClassName('chsp-info-table');
    for (let i = 0; i < infoTables.length; i++) {
        // max height for the chsp-info-table elements is the height of the
        // structure viewer
        (infoTables[i] as HTMLElement).style.maxHeight = `${height}px`;
    }
}

/**
 * Start a download from the user browser, asking it to save the file with the
 * given `filename`. The file will contain the given `content`.
 */
function startDownload(filename: string, content: string): void {
    const a = document.createElement('a');
    a.download = filename;
    a.href = URL.createObjectURL(new Blob([content], { type: 'application/json' }));
    a.style.display = 'none';

    document.body.appendChild(a);
    a.click();
    setTimeout(() => {
        document.body.removeChild(a);
        URL.revokeObjectURL(a.href);
    }, 2000);
}
