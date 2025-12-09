import { Streamlit, RenderData } from 'streamlit-component-lib';
import './chemiscope-streamlit.css';

const ROOT_ID = 'chemiscope-root';
let visualizer: any | null = null;
let indexer: any | null = null;
let visualizerLoaded = false;
let lastSelection: number | null = null;
let lastReportedSelection: number | null = null;
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
        root.style.justifyContent = 'stretch';
    } else if (typeof widthArg === 'number') {
        const viewer = root.querySelector('.chemiscope-streamlit') as HTMLElement;
        if (viewer) {
            viewer.style.width = widthArg + 'px';
        }
        root.style.justifyContent = 'center';
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

function reportSelectionToStreamlit(structureIndex: number | null): void {
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
        let structureIndexToSend: number | null = null;
        
        if (indexes && typeof indexes.structure === 'number') {
            structureIndexToSend = indexes.structure;
        }

        lastSelection = structureIndexToSend;
        reportSelectionToStreamlit(structureIndexToSend);
    };

    if (visualizer.map) {
        originalMapOnselect = visualizer.map.onselect;
        visualizer.map.onselect = (indexes: any) => {
            if (typeof originalMapOnselect === 'function') {
                originalMapOnselect(indexes);
            }
            sendFromIndexes(indexes);
        };
    } else {
        originalMapOnselect = null;
    }

    if (visualizer.structure) {
        originalStructOnselect = visualizer.structure.onselect;
        visualizer.structure.onselect = (indexes: any) => {
            console.log('structure.onselect called with indexes:', indexes);
            if (typeof originalStructOnselect === 'function') {
                originalStructOnselect(indexes);
            }
            sendFromIndexes(indexes);
        };
    } else {
        originalStructOnselect = null;
    }

    if (visualizer.info) {
        originalInfoOnchange = visualizer.info.onchange;
        visualizer.info.onchange = (indexes: any) => {
            console.log('info.onchange called with indexes:', indexes);
            if (typeof originalInfoOnchange === 'function') {
                originalInfoOnchange(indexes);
            }
            sendFromIndexes(indexes);
        };
    } else {
        originalInfoOnchange = null;
    }
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
        displayWarning(ROOT_ID, 'Chemiscope library not loaded. Check script imports.', 0);
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

        const warnings = new Chemiscope.Warnings();
        warnings.addHandler((message: string, timeout: number = 4000) => {
            displayWarning(ROOT_ID, message, timeout);
        });

        toggleLoadingVisible(ROOT_ID, true);

        visualizerClass
            .load(config as any, dataset as any, warnings)
            .then((v: any) => {
                visualizer = v;
                indexer = new Chemiscope.EnvironmentIndexer(dataset.structures, dataset.environments);
                installReverseSyncCallbacks();

                if (selectedIndex !== null) {
                    lastSelection = selectedIndex;
                    applySelectionFromStructure(selectedIndex);
                }
            })
            .catch((err: unknown) => {
                console.error('Error loading visualizer:', err);
                displayWarning(ROOT_ID, 'Error loading visualization: ' + String(err));
            })
            .finally(() => {
                toggleLoadingVisible(ROOT_ID, false);
            });
    } else {
        // 1. Selection Update
        if (selectedIndex !== null && selectedIndex !== lastSelection && visualizer) {
            lastSelection = selectedIndex;
            applySelectionFromStructure(selectedIndex);
        }

        // 2. Settings Update
        if (settings && visualizer) {
            try {
                const newStr = JSON.stringify(settings);
                if (newStr !== JSON.stringify(lastSettings)) {
                    visualizer.applySettings(settings);
                    lastSettings = settings;
                }
            } catch (err) {
                console.error('Error applying settings:', err);
                displayWarning(ROOT_ID, 'Error applying settings: ' + String(err), 4000);
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

function generateHTMLForMode(mode: string, width?: number | 'stretch'): string {
    const widthStyle =
        width && width !== 'stretch' ? `width: ${width}px; margin: 0 auto;` : 'width: 100%;';

    const componentCss = `
  <style>
    #chemiscope-root {
      ${widthStyle}
      height: 100%;
      display: flex;
    }

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

function toggleLoadingVisible(divId: string, visible: boolean = true): void {
    const loader = document.getElementById(`${divId}-loading`);
    if (loader) {
        loader.style.display = visible ? 'block' : 'none';
    }
}

function displayWarning(divId: string, message: any, timeout: number = 4000): void {
    const display = document.getElementById(`${divId}-warning-display`);
    if (!display) {
        console.error('Warning display element not found:', `${divId}-warning-display`);
        return;
    }

    const textMessage = message instanceof Error ? message.toString() : String(message);
    const paragraph = display.getElementsByTagName('p')[0];
    if (paragraph) {
        paragraph.innerText = textMessage;
    }

    display.style.display = 'flex';

    if (timeout > 0) {
        // Automatically remove the warning after the set timeout
        setTimeout(() => {
            display.style.display = 'none';
        }, timeout);
    }
}
