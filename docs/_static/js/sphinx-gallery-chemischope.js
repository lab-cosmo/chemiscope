// function startLoading() {
//     getByID('loading').style.display = 'block';

//     const main = document.getElementsByTagName('main')[0];
//     main.onclick = () => {};
//     main.style.opacity = '0.3';
// }

// function stopLoading() {
//     getByID('loading').style.display = 'none';

//     const main = document.getElementsByTagName('main')[0];
//     main.onclick = null;
//     main.style.opacity = '1';
// }

// function hideLoader() {
//     // Remove loader element
//     var loader = document.getElementById('loader');
//     if (loader) {
//         loader.parentNode.removeChild(loader);
//     }
// }

async function loadChemiscopeSphinxGallery(divId, dataset) {
    // showLoader(divId);
    // setTimeout(() => {}, 1000)
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
    // hideLoader(divId);
};
