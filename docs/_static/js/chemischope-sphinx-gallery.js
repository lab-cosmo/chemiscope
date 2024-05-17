async function loadChemiscopeSphinxGallery(divId, dataset) {
    // Handle warnings
    Chemiscope.addWarningHandler((message) => displayWarning(divId, message));

    // Display loader
    toggleLoaderVisible(divId, true);

    // Load visialisation
    try {
        const config = {
            map: `${divId}-map`,
            info: `${divId}-info`,
            meta: `${divId}-meta`,
            structure: `${divId}-structure`,
        };
        const root = document.getElementById(divId);
        root.innerHTML = generateChemiscopeHTML(config);
        await Chemiscope.DefaultVisualizer.load(config, dataset);
    }

    // Display errors
    catch (error) {
        displayWarning(divId, error)
    }

    // Hide loader
    finally {
        toggleLoaderVisible(divId, false);
    }
}

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

function toggleLoaderVisible(divId, visible = true) {
    const loader = document.getElementById(`${divId}-loading`);
    if (loader) {
        loader.style.display = visible ? "block" : "none";
    }
}

function hideElement(elementId) {
    const element = document.getElementById(elementId);
    if (element) {
        element.style.display = "none";
    } else {
        console.error(`Element ${elementId} is not found`);
    }
}

function displayWarning(divId, message) {
    const display = document.getElementById(`${divId}-warning-display`);
    display.getElementsByTagName('p')[0].innerText = message;
    display.style.display = 'flex';

    // Automatically remove the warning after 4s
    setTimeout(() => {
        display.style.display = 'none';
    }, 4000);
}
