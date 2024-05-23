/**
 * Asynchronously loads the Chemiscope visualization for the dataset
 */
async function loadChemiscopeSphinxGallery(divId, filePath) {
    // Handle warnings
    Chemiscope.addWarningHandler((message) => displayWarning(divId, message));

    // Display loading
    toggleLoadingVisible(divId, true);

    // Load the visualizer
    try {
        const dataset = await fetchDataset(filePath);

        // Setup visualizer config
        const config = {
            map: `${divId}-map`,
            info: `${divId}-info`,
            meta: `${divId}-meta`,
            structure: `${divId}-structure`,
        };

        // Render Chemiscope Visualizer
        const root = document.getElementById(divId);
        root.innerHTML = generateChemiscopeHTML(config);
        await Chemiscope.DefaultVisualizer.load(config, dataset);
    }

    // Display errors
    catch (error) {
        displayWarning(divId, error);
    }

    // Hide loading
    finally {
        toggleLoadingVisible(divId, false);
    }
}

/**
 * Loads the dataset and handled gzipped JSON with NaN values
 */
async function fetchDataset(path) {
    // Load file
    const response = await fetch(path);
    if (!response.ok) {
        throw new Error(`Failed to fetch ${path}: ${response.statusText}`);
    }

    // Get as json
    const buffer = await response.arrayBuffer();
    const decompressedData = pako.inflate(new Uint8Array(buffer), { to: 'string' });
    return parseJsonWithNaN(decompressedData);
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
function generateChemiscopeHTML(config) {
    return `
        <div style="display: flex; flex-wrap: wrap;">
            <div style="flex: 1; padding-right: 10px; position: relative;">
                <div id="${config.meta}"></div>
                <div id="${config.map}"></div>
            </div>
            <div style="flex: 1; padding-left: 10px; position: relative;">
                <div id="${config.structure}" style="width: 100%; height: 100%;"></div>
                <div id="${config.info}"></div>
            </div>
        </div>`;
}

/**
 * Toggles the visibility of the loader spinner
 */
function toggleLoadingVisible(divId, visible = true) {
    const loader = document.getElementById(`${divId}-loading`);
    if (loader) {
        loader.style.display = visible ? "block" : "none";
    }
}

/**
 * Hides the specified element
 */
function hideElement(elementId) {
    const element = document.getElementById(elementId);
    if (element) {
        element.style.display = "none";
    } else {
        console.error(`Element ${elementId} is not found`);
    }
}

/**
 * Displays a warning message in the specified div
 */
function displayWarning(divId, message) {
    const display = document.getElementById(`${divId}-warning-display`);
    display.getElementsByTagName('p')[0].innerText = message;
    display.style.display = 'flex';

    // Automatically remove the warning after 4 seconds
    setTimeout(() => {
        display.style.display = 'none';
    }, 4000);
}
