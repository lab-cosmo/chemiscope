/**
 * Asynchronously loads the Chemiscope visualization for the dataset
 */
async function loadChemiscopeSphinxGallery(divId, dataset) {
    // Handle warnings by displaying them
    Chemiscope.addWarningHandler((message) => displayWarning(divId, message));

    // Display the loading spinner
    toggleLoaderVisible(divId, true);

    try {
        // Configuration for the Chemiscope visualizer
        const config = {
            map: `${divId}-map`,
            info: `${divId}-info`,
            meta: `${divId}-meta`,
            structure: `${divId}-structure`,
        };

        // Get the root element and set its inner HTML to the generated content
        const root = document.getElementById(divId);
        root.innerHTML = generateChemiscopeHTML(config);

        // Load the visualizer with the provided configuration and dataset
        await Chemiscope.DefaultVisualizer.load(config, dataset);
    } catch (error) {
        // Display any errors that occur during loading
        displayWarning(divId, error);
    } finally {
        // Hide the loading spinner
        toggleLoaderVisible(divId, false);
    }
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
function toggleLoaderVisible(divId, visible = true) {
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
