/// Read JSON or gzipped JSON and return the parsed object
function readJSON(path, buffer) {
    let text;
    if (path.endsWith(".gz")) {
        text = pako.inflate(buffer, {to: "string"});
    } else {
        const decoder = new TextDecoder("utf-8");
        text = decoder.decode(buffer);
    }

    // Allow NaN in the JSON file. They are not part of the spec, but
    // Python's json module output them, and they can be useful.
    return JSON.parse(text.replace(/\bNaN\b/g, '"***NaN***"'), function(key, value) {
        return value === "***NaN***" ? NaN : value;
    });
}

let VISUALIZER = undefined;
let DATASET = undefined;
let J2S_PATH = undefined;
let HIDE_ON_DEMAND_STRUCTURES = undefined;

function loadStructureOnDemand(index, structure) {
    /**
     * An example of a loadStructure callback to load structures from an URL on demand
    */
    return JSON.parse($.ajax({
        type: "GET",
        url: structure.data,
        // this is getting deprecated, but the best option for now
        async: false
    }).responseText);
}

function setupChemiscope(json) {
    const config = {
        map:       'chemiscope-map',
        info:      'chemiscope-info',
        meta:      'chemiscope-meta',
        structure: 'chemiscope-structure',
        presets:   json.presets || {},
        j2sPath:   J2S_PATH,
    };

    if (DATASET !== undefined && DATASET.includes("Azaphenacenes.json.gz")) {
        // example of asynchronous structure loading
        config.loadStructure = loadStructureOnDemand;
        HIDE_ON_DEMAND_STRUCTURES.disabled = true;
    } else {
        HIDE_ON_DEMAND_STRUCTURES.disabled = false;
    }

    if (VISUALIZER !== undefined) {
        VISUALIZER.remove();
    }

    Chemiscope.DefaultVisualizer.load(config, json).then((v) => {
        VISUALIZER = v;
        v.structure.positionSettingsModal = (rect) => {
            const structureRect = document.getElementById('chemiscope-structure').getBoundingClientRect();

            return {
                top: structureRect.top,
                left: structureRect.left - rect.width - 25,
            };
        };

        v.map.positionSettingsModal = (rect) => {
            const mapRect = document.getElementById('chemiscope-map').getBoundingClientRect();

            return {
                top: mapRect.top,
                left: mapRect.left + mapRect.width + 25,
            };
        };

        stopLoading();
    }).catch(e => setTimeout(() => {throw e;}));
}

const IGNORED_JSMOL_ERRORS = [
    "IndexSizeError: Failed to execute 'getImageData' on 'CanvasRenderingContext2D': The source width is 0.",
    "IndexSizeError: Index or size is negative or greater than the allowed amount",
];

function displayError(error) {
    if (IGNORED_JSMOL_ERRORS.includes(error.toString())) {
        // Ignores mysterious JSMol resize error we really have no clear way to fix.
        return;
    }

    document.getElementById('loading').style.display = 'none';

    const display = document.getElementById('error-display');
    display.style.display = "block";
    display.getElementsByTagName('p')[0].innerText = error.toString();
    const backtrace = display.getElementsByTagName('details')[0];
    backtrace.getElementsByTagName('p')[0].innerText = error.stack;
}

function displayWarning(message) {
    const display = document.getElementById('warning-display');
    display.style.display = "block";
    display.getElementsByTagName('p')[0].innerText = message;
}

function startLoading() {
    document.getElementById('loading').style.display = 'block';

    const main = document.getElementsByTagName('main')[0];
    main.onclick = () => {};
    main.style.opacity = 0.3;
}

function stopLoading() {
    document.getElementById('loading').style.display = 'none';

    const main = document.getElementsByTagName('main')[0];
    main.onclick = null;
    main.style.opacity = 1;
}

function loadExample(url) {
    startLoading();
    DATASET = url;
    fetch(url)
        .then(response => {
            if (!response.ok) {
                throw Error(`unable to load file at ${url}`)
            } else {
                return response.arrayBuffer();
            }
        })
        .then(buffer => readJSON(url, buffer))
        .then(json => setupChemiscope(json))
        .catch(e => setTimeout(() => {throw e;}));
}

function setupDefaultChemiscope(j2sPath) {
    J2S_PATH = j2sPath;
    Chemiscope.addWarningHandler((message) => displayWarning(message));

    window.onerror = (msg, url, line, col, error) => {
        displayError(error);
    }

    HIDE_ON_DEMAND_STRUCTURES = document.createElement('style');
    HIDE_ON_DEMAND_STRUCTURES.type = 'text/css';
    document.head.appendChild(HIDE_ON_DEMAND_STRUCTURES);
    HIDE_ON_DEMAND_STRUCTURES.sheet.insertRule('.hide-on-demand-structures {display: none}');
    HIDE_ON_DEMAND_STRUCTURES.disabled = false;

    const loadDataset = document.getElementById('load-dataset');
    loadDataset.onchange = () => {
        startLoading();
        const file = loadDataset.files[0];
        DATASET = `user-loaded: ${file.name}`;
        const reader = new FileReader();
        reader.onload = () => {
            const dataset = readJSON(file.name, reader.result);
            setupChemiscope(dataset);
        }
        reader.readAsArrayBuffer(file);
    }

    const saveDataset = document.getElementById('save-dataset');
    const saveDatasetName = document.getElementById('save-dataset-name');
    const includePresets = document.getElementById('save-dataset-presets');
    const includeStructures = document.getElementById('save-dataset-structures');
    saveDataset.onclick = () => {
        const dataset = VISUALIZER.dataset(includeStructures.checked);
        console.log(includePresets.checked);
        if (includePresets.checked) {
            dataset.presets = VISUALIZER.dumpSettings();
        }
        startDownload(saveDatasetName.value, JSON.stringify(dataset));
    }

    const loadPresets = document.getElementById('load-presets');
    loadPresets.onchange = () => {
        const file = loadPresets.files[0];
        const reader = new FileReader();
        reader.onload = () => {
            const presets = readJSON(file.name, reader.result);
            VISUALIZER.applyPresets(presets);
        }
        reader.readAsArrayBuffer(file);
    }

    const savePresets = document.getElementById('save-presets');
    const savePresetsName = document.getElementById('save-presets-name');
    savePresets.onclick = () => {
        startDownload(savePresetsName.value, JSON.stringify(VISUALIZER.dumpSettings()));
    }
}

function startDownload(filename, content) {
    const a = document.createElement('a');
    a.download = filename;
    a.href = URL.createObjectURL(new Blob([content], {type: 'application/json'}));
    a.style.display = 'none';

    document.body.appendChild(a);
    a.click();
    setTimeout(() => {
        document.body.removeChild(a);
        URL.revokeObjectURL(a.href);
    }, 2000);
}
