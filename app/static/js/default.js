function toJson(path, buffer) {
    let text;
    if (path.endsWith(".gz")) {
        text = pako.inflate(buffer, {to: "string"});
    } else {
        const decoder = new TextDecoder("utf-8");
        text = decoder.decode(buffer);
    }

    // Allow NaN in the JSON file. They are not part of the spec, but
    // Python's json module output them, and they can be usefull.
    return JSON.parse(text.replace(/\bNaN\b/g, '"***NaN***"'), function(key, value) {
        return value === "***NaN***" ? NaN : value;
    });
}

let VISUALIZER = undefined;
let DATASET = undefined;
let ORIGINAL_LOAD_STRUCTURE = undefined;

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

function setupChemiscope(dataset) {
    const config = {
        map:       'chemiscope-map',
        info:      'chemiscope-info',
        meta:      'chemiscope-meta',
        structure: 'chemiscope-structure',
        j2sPath:   GLOBAL_J2S_PATH,
    };

    if (DATASET === "Azaphenacenes.json.gz") {
        // example of asynchronous structure loading
        config.loadStructure = loadStructureOnDemand;
    }

    if (VISUALIZER !== undefined) {
        VISUALIZER.remove();
    }

    Chemiscope.DefaultVisualizer.load(config, dataset).then((v) => {
        VISUALIZER = v;
        v.structure.settingsPlacement((rect) => {
            const structureRect = document.getElementById('chemiscope-structure').getBoundingClientRect();

            return {
                top: structureRect.top,
                left: structureRect.left - rect.width - 25,
            };
        });

        v.map.settingsPlacement((rect) => {
            const mapRect = document.getElementById('chemiscope-map').getBoundingClientRect();

            return {
                top: mapRect.top,
                left: mapRect.left + mapRect.width + 25,
            };
        });

        stopLoading();
    }).catch(e => setTimeout(() => {throw e;}));
}

function displayError(error) {
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
        .then(buffer => toJson(url, buffer))
        .then(json => setupChemiscope(json))
        .catch(e => setTimeout(() => {throw e;}));
}

function setupDefaultChemiscope(isStandalone) {
    Chemiscope.addWarningHandler((message) => displayWarning(message));

    window.onerror = (msg, url, line, col, error) => {
        displayError(error);
    }

    const upload = document.getElementById('upload');
    upload.onchange = () => {
        startLoading();
        const name = upload.files[0].name;
        DATASET = `user-loaded: ${name}`;
        const reader = new FileReader();
        reader.onload = () => {
            const json = toJson(name, reader.result);
            setupChemiscope(json);
        }
        reader.readAsArrayBuffer(upload.files[0]);
    }

    if (isStandalone) {
        document.getElementById('examples').style.display = 'none';
        window.GLOBAL_J2S_PATH = 'https://chemapps.stolaf.edu/jmol/jsmol-2019-10-30/j2s/';
    } else {
        window.GLOBAL_J2S_PATH = 'static/js/jsmol/j2s';
    }
}
