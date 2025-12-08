import { Streamlit, RenderData } from 'streamlit-component-lib';

const ROOT_ID = 'chemiscope-root';
let visualizer: any | null = null;
let indexer: any | null = null;
let visualizerLoaded = false;
let lastSelection: number | null = null;
let lastReportedSelection: number | null = null; // CS -> ST
let originalMapOnselect: any | null = null;
let originalStructOnselect: any | null = null;
let originalInfoOnchange: any | null = null;
let lastSettings: any | null = null;

function getChemiscope(): any | null {
    const cs = (window as any).Chemiscope;
    if (!cs) {
        console.error('window.Chemiscope not found. Did chemiscope.min.js load?');
        return null;
    }
    if (!cs.DefaultVisualizer) {
        console.error('Chemiscope.DefaultVisualizer is missing:', cs);
        return null;
    }
    return cs;
}

function getOrCreateRoot(): HTMLDivElement {
    let root = document.getElementById(ROOT_ID) as HTMLDivElement | null;
    if (!root) {
        root = document.createElement('div');
        root.id = ROOT_ID;
        root.style.width = '100%';
        document.body.appendChild(root);
    }
    return root;
}

function applyWidthPolicy(widthArg: string | number, root: HTMLDivElement) {
    if (widthArg === 'stretch') {
        root.style.width = '100%';
    } else if (typeof widthArg === 'number') {
        root.style.width = widthArg + 'px';
        root.style.maxWidth = '100%';
    }
}

function applyHeightPolicy(heightArg: number, root: HTMLDivElement) {
    root.style.height = heightArg + 'px';
    Streamlit.setFrameHeight(heightArg);
}

function applySelectionFromStructure(structureIndex: number): void {
    if (!visualizer || !indexer) {
        return;
    }
    const indexes = indexer.fromStructureAtom('structure', structureIndex);
    if (!indexes) {
        console.warn('No environment for structure index', structureIndex);
        return;
    }

    try {
        visualizer.map.select(indexes);
        if (originalMapOnselect) {
            originalMapOnselect(indexes);
        }
    } catch (err) {
        console.error(err);
    }
}

function reportSelectionToStreamlit(structureIndex: number): void {
    if (structureIndex === lastReportedSelection) {
        return;
    }
    lastReportedSelection = structureIndex;

    try {
        Streamlit.setComponentValue(structureIndex);
    } catch (err) {
        console.error('Error sending selection:', err);
    }
}

function installReverseSyncCallbacks(): void {
    if (!visualizer) return;

    const sendFromIndexes = (indexes: any) => {
        if (!indexes || typeof indexes.structure !== 'number') {
            return;
        }
        lastSelection = indexes.structure;
        reportSelectionToStreamlit(indexes.structure);
    };

    originalMapOnselect = visualizer.map.onselect;
    visualizer.map.onselect = (indexes: any) => {
        if (typeof originalMapOnselect === 'function') {
            originalMapOnselect(indexes);
        }
        sendFromIndexes(indexes);
    };

    originalStructOnselect = visualizer.structure.onselect;
    visualizer.structure.onselect = (indexes: any) => {
        console.log('structure.onselect called with indexes:', indexes);
        if (typeof originalStructOnselect === 'function') {
            originalStructOnselect(indexes);
        }
        sendFromIndexes(indexes);
    };

    // Wrap info.onchange
    originalInfoOnchange = visualizer.info.onchange;
    visualizer.info.onchange = (indexes: any) => {
        console.log('info.onchange called with indexes:', indexes);
        if (typeof originalInfoOnchange === 'function') {
            originalInfoOnchange(indexes);
        }
        sendFromIndexes(indexes);
    };
}

function onRender(event: Event): void {
    const data = (event as CustomEvent<RenderData>).detail;
    const args = data.args as {
        dataset: any;
        height?: number;
        width?: number | string;
        selected_index?: number;
        mode?: string;
        settings?: Map<string, any>;
    };

    const dataset = args.dataset;
    if (!dataset) {
        return;
    }

    const selectedIndex = typeof args.selected_index === 'number' ? args.selected_index : null;
    const mode = typeof args.mode === 'string' ? args.mode : 'default';
    const settings = args.settings || {};
    const widthArg = args.width ?? 'stretch';
    const heightArg = typeof args.height === 'number' ? args.height : 550;

    const Chemiscope = getChemiscope();
    if (!Chemiscope) {
        return;
    }

    // First render
    if (!visualizerLoaded) {
        visualizerLoaded = true;

        const root = getOrCreateRoot();
        root.innerHTML = generateHTMLForMode(mode);

        applyWidthPolicy(widthArg, root);
        applyHeightPolicy(heightArg, root);

        const config = {
            map: 'chemiscope-map',
            info: 'chemiscope-info',
            meta: 'chemiscope-meta',
            structure: 'chemiscope-structure',
        };
        const visualizerClass = getVisualizerForMode(mode, Chemiscope);

        if (Object.keys(settings).length > 0) {
            try {
                dataset.settings = Object.assign({}, dataset.settings, settings);
            } catch {
                console.warn('Could not attach settings to dataset');
            }
        }

        visualizerClass
            .load(config as any, dataset as any)
            .then((v: any) => {
                visualizer = v;
                const ds = v.dataset(true);
                indexer = new Chemiscope.EnvironmentIndexer(ds.structures, ds.environments);
                installReverseSyncCallbacks();

                if (selectedIndex !== null) {
                    lastSelection = selectedIndex;
                    applySelectionFromStructure(selectedIndex);
                }
            })
            .catch((err: unknown) => console.error('Error loading visualizer:', err));
    } else {
        if (selectedIndex !== null && selectedIndex !== lastSelection && visualizer) {
            lastSelection = selectedIndex;
            applySelectionFromStructure(selectedIndex);
        }

        if (settings && visualizer) {
            try {
                const newStr = JSON.stringify(settings);
                if (newStr !== JSON.stringify(lastSettings)) {
                    visualizer.applySettings(settings);
                    lastSettings = settings;
                }
            } catch (err) {
                console.error('Error applying settings:', err);
            }
        }

        const root = getOrCreateRoot();
        applyWidthPolicy(widthArg, root);
        applyHeightPolicy(heightArg, root);
    }
}

Streamlit.events.addEventListener(Streamlit.RENDER_EVENT, onRender);
Streamlit.setComponentReady();

function getVisualizerForMode(mode: string, Chemiscope: any): any {
    switch (mode) {
        case 'structure':
            return Chemiscope.StructureVisualizer;
        case 'map':
            return Chemiscope.MapVisualizer;
        default:
            return Chemiscope.DefaultVisualizer;
    }
}

function generateHTMLForMode(mode: string): string {
    const componentCss = `
  <style>
    .chemiscope-streamlit { width: 100%; height: 100%; }
    .visualizer-container { display: flex; flex-direction: row; width: 100%; height: 100%; }
    .visualizer-column, .visualizer-column-right {
      flex: 1 1 0; width: 50%; height: 100%; display: flex; flex-direction: column;
    }
    .visualizer-item { flex: 1; width: 100%; height: 100%; }
    .visualizer-info { flex-shrink: 0; }
    .visualizer-map-mode .visualizer-column,
    .visualizer-structure-mode .visualizer-column { width: 100%; flex: 1 1 100%; }
  </style>`;

    switch (mode) {
        case 'structure':
            return `
      ${componentCss}
      <div class="chemiscope-streamlit">
        <div class="visualizer-container visualizer-structure-mode">
          <div class="visualizer-column">
            <div id="chemiscope-meta"></div>
            <div id="chemiscope-structure" class="visualizer-item"></div>
            <div id="chemiscope-info" class="visualizer-info"></div>
          </div>
        </div>
      </div>`;

        case 'map':
            return `
      ${componentCss}
      <div class="chemiscope-streamlit">
        <div class="visualizer-container visualizer-map-mode">
          <div class="visualizer-column">
            <div id="chemiscope-meta"></div>
            <div id="chemiscope-map" class="visualizer-item"></div>
            <div id="chemiscope-info" class="visualizer-info"></div>
          </div>
        </div>
      </div>`;

        default:
            return `
      ${componentCss}
      <div class="chemiscope-streamlit">
        <div class="visualizer-container">
          <div class="visualizer-column-right">
            <div id="chemiscope-meta"></div>
            <div id="chemiscope-map" class="visualizer-item"></div>
          </div>
          <div class="visualizer-column">
            <div id="chemiscope-structure" class="visualizer-item"></div>
            <div id="chemiscope-info" class="visualizer-info"></div>
          </div>
        </div>
      </div>`;
    }
}
