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
    selected_index?: number | undefined;
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
        select: (indexes: any) => void;
        onselect: ((indexes: any) => void) | null;
        activeChanged?: ((guid: any, indexes: any) => void) | null;
    };
    info?: {
        onchange: ((indexes: any) => void) | null;
    };
    select: (indexes: any) => void;
    applySettings: (settings: Record<string, any>) => void;
    saveSettings: () => Record<string, any>;
    onSettingChange: (callback: (keys: string[], value: unknown) => void) => void;
}

enum StreamlitValue {
    SELECTION = 'selected_id',
    SETTINGS = 'settings',
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

        // Track selections
        currentSelection: null as number | null,

        // Track settings
        currentSettings: null as string | null,
        
        // Track source of the selection events to avoid loops
        selectionFromPython: false,
        settingsFromPython: false, 

        // Original callbacks
        originalMapOnselect: null as any,
        originalStructOnselect: null as any,
        originalStructActiveChanged: null as any,
        originalSelect: null as any,
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

        const mode = (typeof args.mode === 'string' ? args.mode : 'default') as ChemiscopeMode;
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
                args.settings || {},
                args.selected_index,
                widthArg,
                heightArg
            );
        } else {
            this.handleUpdate(args, widthArg, heightArg);
        }
    }

    private handleFirstRender(
        Chemiscope: any,
        dataset: Record<string, any>,
        mode: ChemiscopeMode,
        settings: Record<string, any>,
        selectedIndex: number | undefined,
        widthArg: string | number,
        heightArg: number
    ): void {
        this.state.loaded = true;

        const root = getOrCreateRoot();
        root.innerHTML = generateHTMLForMode(mode);

        applyWidthPolicy(widthArg, root);
        applyHeightPolicy(heightArg, root);

        // Store initial received values
        this.state.currentSelection = selectedIndex !== undefined ? selectedIndex : null;
        this.state.currentSettings = JSON.stringify(settings);

        this.initializeVisualizer(Chemiscope, dataset, settings, selectedIndex);
    }

    private handleUpdate(args: ChemiscopeArgs, widthArg: string | number, heightArg: number): void {
        const root = getOrCreateRoot();
        applyWidthPolicy(widthArg, root);
        applyHeightPolicy(heightArg, root);

        // Handle traitlet-style updates
        this.handleSettingsUpdate(args.settings);
        this.handleSelectionUpdate(args.selected_index);
    }

    private handleSelectionUpdate(selectedIndex: number | null | undefined): void {
        console.log('Handling selection update:', selectedIndex);
        if (selectedIndex === undefined) {
            return;
        }
        
        // Only process if selection actually changed
        if (selectedIndex !== this.state.currentSelection) {
            console.log("update selection");
            
            this.state.selectionFromPython = true;
            try {
                this.state.currentSelection = selectedIndex;
                this.applySelection(selectedIndex);
            } finally {
                this.state.selectionFromPython = false;
            }
        }
    }

    private handleSettingsUpdate(settings: Record<string, any> | undefined): void {
        console.log('Handling settings update')
        if (!settings || !this.state.visualizer) {
            return;
        }

        const settingsStr = JSON.stringify(settings);

        // Only process if settings actually changed
        if (settingsStr !== this.state.currentSettings) {
            console.log("update settings");            

            this.state.settingsFromPython = true;
            try {
                this.state.currentSettings = settingsStr;
                this.state.visualizer.applySettings(settings);
            } finally {
                this.state.settingsFromPython = false;
            }
        }
    }

    private applySelection(selectedIndex: number | null | undefined): void {
        console.log('Applying selection:', selectedIndex, 'from python', this.state.selectionFromPython);
        const { visualizer, indexer } = this.state;
        if (!visualizer || !indexer) {
            return;
        }

        if (selectedIndex === null) {
            // Handle deselection if needed
            return;
        }

        const indexes = indexer.fromStructureAtom('structure', selectedIndex);
        if (!indexes) {
            console.warn('No environment for structure index', selectedIndex);
            return;
        }
        
        console.log("calling originalselect, ", indexes);
        this.state.originalSelect?.(indexes);
    }

    private sendSelectionToStreamlit(indexes: any): void {
        console.log('Sending selection to Streamlit:', indexes, this.state.selectionFromPython);
        if (this.state.selectionFromPython) {
            return;
        }

        let structureIdToSend: number | null = null;
        if (indexes && typeof indexes.structure === 'number') {
            structureIdToSend = indexes.structure;
        }

        // Only send if changed
        if (structureIdToSend !== this.state.currentSelection) {
            this.state.currentSelection = structureIdToSend;
            
            // Get current settings
            const currentSettings = this.state.visualizer?.saveSettings() || {};
            this.state.currentSettings = JSON.stringify(currentSettings);

            Streamlit.setComponentValue({
                [StreamlitValue.SETTINGS]: currentSettings,
                [StreamlitValue.SELECTION]: structureIdToSend,
            });
        }
    }

    private sendSettingsToStreamlit(settings: Record<string, any>): void {
        console.log('Sending settings to Streamlit:', settings, this.state.settingsFromPython);        
        if (this.state.settingsFromPython) {
            return;
        }

        const settingsStr = JSON.stringify(settings);

        // Only send if changed
        if (settingsStr !== this.state.currentSettings) {
            this.state.currentSettings = settingsStr;
            
            Streamlit.setComponentValue({
                [StreamlitValue.SETTINGS]: settings,
                [StreamlitValue.SELECTION]: this.state.currentSelection,
            });
        }
    }

    private installReverseSyncCallbacks(): void {
        const visualizer = this.state.visualizer;
        if (!visualizer) {
            return;
        }

        // Map onselect
        if (visualizer.map) {
            const originalMapOnselect = visualizer.map.onselect?.bind(visualizer.map);

            this.state.originalMapOnselect = originalMapOnselect;
            visualizer.map.onselect = (indexes: any) => {
                console.log("map.onselect");
                originalMapOnselect?.(indexes);
                this.sendSelectionToStreamlit(indexes);
            };
        }

        // Structure onselect - selection from structure viewer to Streamlit
        if (visualizer.structure) {
            const originalStructOnselect = visualizer.structure.onselect?.bind(
                visualizer.structure
            );
            this.state.originalStructOnselect = originalStructOnselect;
            visualizer.structure.onselect = (indexes: any) => {
                console.log("structure.onselect");
                originalStructOnselect?.(indexes);
                this.sendSelectionToStreamlit(indexes);
            };

            if (typeof visualizer.structure.activeChanged === 'function') {
                const originalActiveChanged = visualizer.structure.activeChanged.bind(
                    visualizer.structure
                );
                this.state.originalStructActiveChanged = originalActiveChanged;

                visualizer.structure.activeChanged = (guid: any, indexes: any) => {
                    console.log("structure.activeChanged");
                    originalActiveChanged?.(guid, indexes);
                    this.sendSelectionToStreamlit(indexes);
                };
            }
        }

        // Info onchange
        if (visualizer.info) {
            const originalSelect = visualizer.select.bind(visualizer);
            this.state.originalSelect = originalSelect;

            visualizer.select = (indexes: any) => {
                console.log("visualizer.select");
                 
                originalSelect(indexes);
                this.sendSelectionToStreamlit(indexes);
            };
        }

        // Settings change - from visualizer to Streamlit
        visualizer.onSettingChange((_keys: string[], _value: unknown) => {
            console.log("onsettingschange");
            const currentSettings = visualizer.saveSettings();
            this.sendSettingsToStreamlit(currentSettings);
        });
    }

    private initializeVisualizer(
        Chemiscope: any,
        dataset: Record<string, any>,
        settings: Record<string, any>,
        selectedIndex: number | null | undefined
    ): void {
        const mode = dataset.metadata?.mode || 'default';
        const visualizerClass = getVisualizerClass(mode, Chemiscope);

        // Merge initial settings
        try {
            dataset.settings = Object.assign({}, dataset.settings, settings);
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

                // Apply initial selection if provided
                if (selectedIndex !== null && selectedIndex !== undefined) {
                    this.state.selectionFromPython = true;
                    try {
                        this.state.currentSelection = selectedIndex;
                        this.applySelection(selectedIndex);                    
                    } finally {
                        this.state.selectionFromPython = false;
                    }
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
}