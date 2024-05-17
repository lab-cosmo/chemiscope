async function loadChemiscopeSphinxGallery(divId, dataset) {
    showLoader(divId);
    const config = {
        map: `${divId}-map`,
        info: `${divId}-info`,
        meta: `${divId}-meta`,
        structure: `${divId}-structure`,
    };
    const root = document.getElementById(divId);
    root.innerHTML = generateChemiscopeHTML(config);

    try {
        await Chemiscope.DefaultVisualizer.load(config, dataset);
    } catch (error) {
        console.error(error);
    } finally {
        hideLoader(divId);
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

function toggleLoader(divId, show) {
    const loader = document.getElementById(`${divId}-loading`);
    if (loader) {
        loader.style.display = show ? "block" : "none";
    }
}
