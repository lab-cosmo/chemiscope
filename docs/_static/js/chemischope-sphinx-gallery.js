async function loadChemiscopeSphinxGallery(divId, dataset) {
    showLoader();
    const root = document.getElementById(divId);
    root.innerHTML = `
        <div style="display: flex; flex-wrap: wrap;">
            <div style="flex: 1; padding-right: 10px; position: relative;">
                <div id="chemiscope-meta"></div>
                <div id="chemiscope-map"></div>
            </div>
            <div style="flex: 1; padding-left: 10px; position: relative;">
                <div id="chemiscope-structure" style="width: 100%; height: 100%;"></div>
                <div id="chemiscope-info"></div>
            </div>
        </div>`;
    const config = {
        map: 'chemiscope-map',
        info: 'chemiscope-info',
        meta: 'chemiscope-meta',
        structure: 'chemiscope-structure',
    };
    const visualizer = await Chemiscope.DefaultVisualizer.load(config, dataset);
    hideLoader();
};

function showLoader() {
    const loader = document.getElementById('loading');
    loader.style.display = 'block';
}

function hideLoader() {
    const loader = document.getElementById('loading');
    loader.style.display = 'none';
}
