// load bootstrap: this needs to come first so that CSS files are loaded in the
// right order
import 'bootstrap';
import 'bootstrap/dist/css/bootstrap.min.css';

import { getByID, logger } from '../src/utils';
import { Dataset, Structure } from '../src/dataset';
import { version, DefaultVisualizer, Settings } from '../src/index';

import { inflate } from 'pako';

// load CSS for the app
import './app.css';

interface Configuration {
    /// optional callback to load the structures on demand.
    loadStructure?: (index: number, structure: unknown) => Structure;
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
        logger.addHandler('warn', (message) => displayWarning(message as string));

        // when the window is resized, change the size available to the info
        // widget
        window.addEventListener('resize', updateInfoWidgetHeight);

        // setup the main HTML
        const root = getByID(id);
        root.innerHTML = `<div class="container-fluid">
            <div class="row">
                <div class="col-md-6" style="padding: 0">
                    <div class="ratio ratio-1x1">
                        <div id="chemiscope-meta" style="z-index: 10"></div>
                        <div id="chemiscope-map" style="position: absolute"></div>
                    </div>
                </div>

                <div class="col-md-6" style="padding: 0">
                    <div class="ratio ratio-5x7">
                        <div>
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
        </div>`;

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
        let config: Partial<Configuration> = {
            loadStructure: undefined,
        };
        if (example == 'Azaphenacenes') {
            // example of asynchronous structure loading
            config.loadStructure = (_: number, structure: any) => {
                return JSON.parse(
                    $.ajax({
                        type: 'GET',
                        url: `examples/${structure.data}`,
                        // this is getting deprecated, but the best option until
                        // we change loadStructure to be an async callback
                        async: false,
                    }).responseText
                );
            };
        }

        await this.fetchAndLoad(`/examples/${example}.json.gz`, config);
    }

    /**
     * Fetch the JSON file at the given URL, an load it as a chemiscope dataset.
     *
     * Optionally specify the configuration in `config`
     */
    public async fetchAndLoad(url: string, config: Partial<Configuration> = {}): Promise<void> {
        startLoading();
        this.dataset = url;

        const response = await fetch(this.dataset);
        if (!response.ok) {
            throw Error(`unable to load file at '${this.dataset}'`);
        }
        const buffer = await response.arrayBuffer();

        if (buffer.byteLength === 0) {
            // Until https://github.com/gnuns/allOrigins/issues/70 is resolved,
            // we can not check if the URL does not exists, so let's assume an
            // error if the buffer is empty
            if (this.dataset.startsWith('https://api.allorigins.win/raw?url=')) {
                let originalUrl = decodeURIComponent(this.dataset.substr(35));
                throw Error(`unable to load file at '${originalUrl}'`);
            }
        }

        const dataset = readJSON(this.dataset, buffer);
        await this.load(config as Configuration, dataset);
    }

    /**
     * Load the given `dataset` with the specified `configuration`.
     */

    public async load(configuration: Configuration, dataset: Dataset): Promise<void> {
        // hide any error coming from the previous dataset loading
        const errors = getByID('error-display');
        errors.style.display = 'none';

        const config = {
            map: 'chemiscope-map',
            info: 'chemiscope-info',
            meta: 'chemiscope-meta',
            structure: 'chemiscope-structure',
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
            let left = structureRect.left - rect.width + 25;
            if (left < 25) {
                left = 25;
            }
            return {
                top: structureRect.top,
                left: left,
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
                this.load({}, dataset);
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
    display.getElementsByTagName('p')[0].innerText = message;
    display.style.display = 'block';

    // automatically remove the warning after 4s
    setTimeout(() => {
        display.style.display = 'none';
    }, 4000);
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
