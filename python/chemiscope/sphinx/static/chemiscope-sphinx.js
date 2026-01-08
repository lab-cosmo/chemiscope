/* eslint-disable */

/**
 * Enum of the modes
 */
const VISUALISER_MODE = {
    DEFAULT: 'default',
    STRUCTURE: 'structure',
    MAP: 'map',
};

/**
 * Asynchronously loads the Chemiscope visualization for the dataset
 */
async function loadChemiscopeSphinx(
    divId,
    filePath,
    visualizerMode = VISUALISER_MODE.DEFAULT,
    timeout = 4000
) {
    // Display loading
    toggleLoadingVisible(divId, true);

    const warnings = new Chemiscope.Warnings();
    warnings.defaultTimeout = timeout; // defaults to 4s visibility
    warnings.addHandler((message, timeout) => {
        displayWarning(divId, message, timeout);
    });

    // Load the visualizer
    try {
        const dataset = await fetchDataset(filePath);

        // Setup visualizer config
        let config = {
            map: `${divId}-map`,
            info: `${divId}-info`,
            meta: `${divId}-meta`,
            structure: `${divId}-structure`,
        };

        const has_external_structures =
            dataset.structures && dataset.structures.length > 0 && dataset.structures[0].data;
        if (visualizerMode !== VISUALISER_MODE.MAP && has_external_structures) {
            // base href for external structures
            const baseHref = filePath + '-ext/';

            config.loadStructure = async (_, structure) => {
                const url = baseHref + String(structure.data).replace(/^\/+/, '');

                return await fetchDataset(url);
            };
        }

        // Prepare html for the visualizer
        const root = document.getElementById(divId);
        root.innerHTML = generateChemiscopeHTML(config, visualizerMode);

        // Load widget
        const visualiser = getVisualizer(visualizerMode);
        await visualiser.load(config, dataset, warnings);
    } catch (error) {
        // Display errors
        console.error(error);
        displayWarning(divId, error);
    } finally {
        // Hide loading
        toggleLoadingVisible(divId, false);
    }
}

/**
 * Returns the appropriate Chemiscope visualizer based on the given mode
 */
function getVisualizer(mode) {
    switch (mode) {
        case VISUALISER_MODE.STRUCTURE:
            return Chemiscope.StructureVisualizer;
        case VISUALISER_MODE.MAP:
            return Chemiscope.MapVisualizer;
        default:
            return Chemiscope.DefaultVisualizer;
    }
}

/**
 * Loads the dataset and handled gzipped JSON with NaN values
 */
async function fetchDataset(filePath) {
    // Get path to the file
    const baseUrl = `${filePath}`;

    // Load file
    const response = await fetch(baseUrl);
    if (!response.ok) {
        throw new Error(`Failed to fetch ${filePath}: ${response.statusText}`);
    }
    const buffer = await response.arrayBuffer();
    const magic = new Uint8Array(buffer.slice(0, 2));

    let text;
    // '1f 8b' is the magic constant starting gzip files
    if (magic[0] == 0x1f && magic[1] == 0x8b) {
        text = pako.inflate(new Uint8Array(buffer), { to: 'string' });
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
function parseJsonWithNaN(text) {
    return JSON.parse(text.replace(/\bNaN\b/g, '"***NaN***"'), (_key, value) => {
        return value === '***NaN***' ? NaN : value;
    });
}

/**
 * Generates the HTML content for the Chemiscope visualizer
 */
function generateChemiscopeHTML(config, visualizerMode) {
    switch (visualizerMode) {
        case VISUALISER_MODE.DEFAULT:
            return `
                <div class="chemiscope-sphinx">
                    <div class="visualizer-container">
                        <div class="visualizer-column-right">
                            <div id="${config.meta}"></div>
                            <div id="${config.map}"></div>
                        </div>
                        <div class="visualizer-column">
                            <div id="${config.structure}" class="visualizer-item"></div>
                            <div id="${config.info}" class="visualizer-info"></div>
                        </div>
                    </div>
                </div>`;
        case VISUALISER_MODE.STRUCTURE:
            return `
                <div class="chemiscope-sphinx">
                    <div class="visualizer-container visualizer-structure-mode">
                        <div class="visualizer-column">
                            <div id="${config.meta}"></div>
                            <div id="${config.structure}" class="visualizer-item"></div>
                            <div id="${config.info}" class="visualizer-info"></div>
                        </div>
                    </div>
                </div>`;
        case VISUALISER_MODE.MAP:
            return `
                <div class="chemiscope-sphinx">
                <div class="visualizer-container visualizer-map-mode">
                    <div class="visualizer-column">
                        <div id="${config.meta}"></div>
                        <div id="${config.map}" class="visualizer-item"></div>
                        <div id="${config.info}" class="visualizer-info"></div>
                    </div>
                </div>
            </div>`;
    }
}

/**
 * Toggles the visibility of the loader spinner
 */
function toggleLoadingVisible(divId, visible = true) {
    const loader = document.getElementById(`${divId}-loading`);
    if (loader) {
        loader.style.display = visible ? 'block' : 'none';
    }
}

/**
 * Hides the specified element
 */
function hideElement(elementId) {
    const element = document.getElementById(elementId);
    if (element) {
        element.style.display = 'none';
    } else {
        console.error(`Element ${elementId} is not found`);
    }
}

/**
 * Displays a warning message in the specified div
 */
function displayWarning(divId, message, timeout) {
    if (timeout < 0) {
        return;
    }
    const display = document.getElementById(`${divId}-warning-display`);
    display.getElementsByTagName('p')[0].innerText = message;
    display.style.display = 'flex';

    if (timeout > 0) {
        // Automatically remove the warning after the set timeout seconds
        setTimeout(() => {
            display.style.display = 'none';
        }, timeout);
    }
}
