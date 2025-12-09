import { RenderData, Streamlit } from 'streamlit-component-lib';
import {
    CONFIG,
    applyHeightPolicy,
    applyWidthPolicy,
    displayWarning,
    generateHTMLForMode,
    getOrCreateRoot,
    toggleLoadingVisible,
} from './dom-utils';

type ChemiscopeMode = 'default' | 'structure' | 'map';

interface ChemiscopeArgs {
    dataset: Record<string, any>;
    height?: number;
    width?: number | string;
    selected_index?: number;
    mode?: ChemiscopeMode;
    settings?: Record<string, any>;
}

interface ChemiscopeWindow extends Window {
    Chemiscope?: any;
}

interface ChemiscopeVisualizer {
    map?: {
        select: (indexes: any) => void;
        onselect: ((indexes: any) => void) | null;
    };
    structure?: {
        onselect: ((indexes: any) => void) | null;
    };
    info?: {
        onchange: ((indexes: any) => void) | null;
    };
    applySettings: (settings: Record<string, any>) => void;
}

function getChemiscope(): any | null {
    const cs = (window as unknown as ChemiscopeWindow).Chemiscope;
    if (!cs) {
        console.error('window.Chemiscope not found. Did chemiscope.min.js load?');
        return null;
    }
    return cs;
}

function getVisualizerClass(mode: string, Chemiscope: any): any {
    switch (mode) {
        case 'structure':
            return Chemiscope.StructureVisualizer;
        case 'map':
            return Chemiscope.MapVisualizer;
        default:
            return Chemiscope.DefaultVisualizer;
    }
}

export class ChemiscopeComponent {
    private state = {
        visualizer: null as ChemiscopeVisualizer | null,
        indexer: null as any | null,
        loaded: false,
        lastSelection: null as number | null,
        lastReportedSelection: null as number | null,
        lastSettings: null as string | null,
        originalMapOnselect: null as any,
        originalStructOnselect: null as any,
        originalInfoOnchange: null as any,
    };

    constructor() {
        Streamlit.events.addEventListener(Streamlit.RENDER_EVENT, (event: Event) => {
            this.onRender((event as CustomEvent<RenderData>).detail);
        });
        Streamlit.setComponentReady();
    }

    private onRender(data: RenderData): void {
        const args = data.args as ChemiscopeArgs;

        const dataset = args.dataset;
        if (!dataset) {
            return;
        }

        const selectedIndex = typeof args.selected_index === 'number' ? args.selected_index : null;
        const mode: ChemiscopeMode = (
            typeof args.mode === 'string' ? args.mode : 'default'
        ) as ChemiscopeMode;
        const settings = args.settings || {};
        const widthArg = args.width ?? 'stretch';
        const heightArg = typeof args.height === 'number' ? args.height : 550;

        const Chemiscope = getChemiscope();
        if (!Chemiscope) {
            displayWarning('Chemiscope library not loaded. Check script imports.', 0);
            return;
        }

        if (!this.state.loaded) {
            this.handleFirstRender(
                Chemiscope,
                dataset,
                mode,
                settings,
                selectedIndex,
                widthArg,
                heightArg
            );
        } else {
            this.handleUpdate(selectedIndex, settings, widthArg, heightArg);
        }
    }

    private handleFirstRender(
        Chemiscope: any,
        dataset: Record<string, any>,
        mode: ChemiscopeMode,
        settings: Record<string, any>,
        selectedIndex: number | null,
        widthArg: string | number,
        heightArg: number
    ): void {
        this.state.loaded = true;

        const root = getOrCreateRoot();
        root.innerHTML = generateHTMLForMode(mode);

        applyWidthPolicy(widthArg, root);
        applyHeightPolicy(heightArg, root);

        this.initializeVisualizer(Chemiscope, dataset, settings, selectedIndex);
    }

    private handleUpdate(
        selectedIndex: number | null,
        settings: Record<string, any>,
        widthArg: string | number,
        heightArg: number
    ): void {
        // Selection update
        if (selectedIndex !== null && selectedIndex !== this.state.lastSelection) {
            this.state.lastSelection = selectedIndex;
            this.applySelectionFromStructure(selectedIndex);
        }

        // Settings update
        this.handleSettingsUpdate(settings);

        const root = getOrCreateRoot();
        applyWidthPolicy(widthArg, root);
        applyHeightPolicy(heightArg, root);
    }

    private applySelectionFromStructure(structureIndex: number): void {
        const { visualizer, indexer } = this.state;
        if (!visualizer || !indexer) {
            return;
        }

        const indexes = indexer.fromStructureAtom('structure', structureIndex);
        if (!indexes) {
            console.warn('No environment for structure index', structureIndex);
            return;
        }

        try {
            visualizer.map?.select(indexes);
            this.state.originalMapOnselect?.(indexes);
        } catch (err) {
            console.error('Error applying selection:', err);
        }
    }

    private sendSelectionToStreamlit(indexes: any): void {
        let structureIdToSend: number | null = null;

        if (indexes && typeof indexes.structure === 'number') {
            structureIdToSend = indexes.structure;
        }

        if (structureIdToSend === this.state.lastSelection) {
            return;
        }
        this.state.lastSelection = structureIdToSend;

        if (structureIdToSend === this.state.lastReportedSelection) {
            return;
        }
        this.state.lastReportedSelection = structureIdToSend;

        try {
            Streamlit.setComponentValue(structureIdToSend);
        } catch (err) {
            console.error('Error sending selection:', err);
        }
    }

    private installReverseSyncCallbacks(): void {
        const visualizer = this.state.visualizer;
        if (!visualizer) {
            return;
        }

        // Map onselect
        if (visualizer.map) {
            this.state.originalMapOnselect = visualizer.map.onselect;
            visualizer.map.onselect = (indexes: any) => {
                this.state.originalMapOnselect?.(indexes);
                this.sendSelectionToStreamlit(indexes);
            };
        } else {
            this.state.originalMapOnselect = null;
        }

        // Structure onselect
        if (visualizer.structure) {
            this.state.originalStructOnselect = visualizer.structure.onselect;
            visualizer.structure.onselect = (indexes: any) => {
                this.state.originalStructOnselect?.(indexes);
                this.sendSelectionToStreamlit(indexes);
            };
        } else {
            this.state.originalStructOnselect = null;
        }

        // Info onchange
        if (visualizer.info) {
            this.state.originalInfoOnchange = visualizer.info.onchange;
            visualizer.info.onchange = (indexes: any) => {
                console.log('info.onchange called with indexes:', indexes);
                this.state.originalInfoOnchange?.(indexes);
                this.sendSelectionToStreamlit(indexes);
            };
        } else {
            this.state.originalInfoOnchange = null;
        }
    }

    private initializeVisualizer(
        Chemiscope: any,
        dataset: Record<string, any>,
        settings: Record<string, any>,
        selectedIndex: number | null
    ): void {
        const mode = dataset.metadata?.mode || 'default';
        const visualizerClass = getVisualizerClass(mode, Chemiscope);

        // Merge user settings into dataset settings
        try {
            dataset.settings = Object.assign({}, dataset.settings, settings);
            this.state.lastSettings = JSON.stringify(settings);
        } catch (e) {
            console.warn('Could not attach settings to dataset:', e);
        }

        const warnings = new Chemiscope.Warnings();
        warnings.addHandler((message: string, timeout: number = 4000) => {
            displayWarning(message, timeout);
        });

        toggleLoadingVisible(true);

        visualizerClass
            .load(CONFIG, dataset, warnings)
            .then((v: ChemiscopeVisualizer) => {
                this.state.visualizer = v;
                this.state.indexer = new Chemiscope.EnvironmentIndexer(
                    dataset.structures,
                    dataset.environments
                );
                this.installReverseSyncCallbacks();

                if (selectedIndex !== null) {
                    this.state.lastSelection = selectedIndex;
                    this.applySelectionFromStructure(selectedIndex);
                }
            })
            .catch((err: unknown) => {
                console.error('Error loading visualizer:', err);
                displayWarning('Error loading visualization: ' + String(err));
            })
            .finally(() => {
                toggleLoadingVisible(false);
            });
    }

    private handleSettingsUpdate(settings: Record<string, any>): void {
        if (!this.state.visualizer) {
            return;
        }

        try {
            const newStr = JSON.stringify(settings);
            if (newStr !== this.state.lastSettings) {
                this.state.visualizer.applySettings(settings);
                this.state.lastSettings = newStr;
            }
        } catch (err) {
            console.error('Error applying settings:', err);
            displayWarning('Error applying settings: ' + String(err), 4000);
        }
    }
}
