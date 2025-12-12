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
        lastSentSelection: null as number | null,
        lastReceivedSelection: null as number | null,

        // Track settings
        lastSentSettings: null as string | null,
        lastReceivedSettings: null as string | null,

        // Prevent update loops
        isProcessingSelection: false,
        isProcessingSettings: false,

        // Original callbacks
        originalMapOnselect: null as any,
        originalStructOnselect: null as any,
        originalStructActiveChanged: null as any,
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
        this.state.lastReceivedSelection = selectedIndex !== undefined ? selectedIndex : null;
        this.state.lastReceivedSettings = JSON.stringify(settings);

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
        if (selectedIndex === undefined) {
            return;
        }

        // Only process if selection actually changed
        if (selectedIndex !== this.state.lastReceivedSelection) {
            this.state.lastReceivedSelection = selectedIndex;

            // Skip if we're already processing a selection
            if (!this.state.isProcessingSelection) {
                this.applyExternalSelection(selectedIndex);
            }
        }
    }

    private handleSettingsUpdate(settings: Record<string, any> | undefined): void {
        if (!settings || !this.state.visualizer) {
            return;
        }

        const settingsStr = JSON.stringify(settings);

        // Only process if settings actually changed
        if (settingsStr !== this.state.lastReceivedSettings) {
            this.state.lastReceivedSettings = settingsStr;

            // Skip if we're already processing settings
            if (!this.state.isProcessingSettings) {
                this.state.isProcessingSettings = true;
                try {
                    this.state.visualizer.applySettings(settings);
                } catch (err) {
                    console.error('Error applying settings:', err);
                } finally {
                    this.state.isProcessingSettings = false;
                }
            }
        }
    }

    private applyExternalSelection(selectedIndex: number | null | undefined): void {
        const { visualizer, indexer } = this.state;
        if (!visualizer || !indexer) {
            return;
        }

        this.state.isProcessingSelection = true;

        try {
            if (selectedIndex === null) {
                // Handle deselection if needed
                return;
            }

            const indexes = indexer.fromStructureAtom('structure', selectedIndex);
            if (!indexes) {
                console.warn('No environment for structure index', selectedIndex);
                return;
            }

            if (visualizer.map) {
                visualizer.map.select(indexes);
            }
            if (visualizer.structure && typeof visualizer.structure.select === 'function') {
                visualizer.structure.select(indexes);
            }
            if (visualizer.info && typeof visualizer.info.onchange === 'function') {
                visualizer.info.onchange(indexes);
            }
        } finally {
            setTimeout(() => {
                this.state.isProcessingSelection = false;
            }, 100);
        }
    }

    private sendSelectionToStreamlit(indexes: any): void {
        if (this.state.isProcessingSelection) {
            return;
        }

        let structureIdToSend: number | null = null;
        if (indexes && typeof indexes.structure === 'number') {
            structureIdToSend = indexes.structure;
        }

        // Only send if changed
        if (structureIdToSend !== this.state.lastSentSelection) {
            this.state.lastSentSelection = structureIdToSend;
            this.state.isProcessingSelection = true;

            // Get current settings
            const currentSettings = this.state.visualizer?.saveSettings() || {};

            Streamlit.setComponentValue({
                [StreamlitValue.SETTINGS]: currentSettings,
                [StreamlitValue.SELECTION]: structureIdToSend,
            });

            setTimeout(() => {
                this.state.isProcessingSelection = false;
            }, 100);
        }
    }

    private sendSettingsToStreamlit(settings: Record<string, any>): void {
        if (this.state.isProcessingSettings) {
            return;
        }

        const settingsStr = JSON.stringify(settings);

        // Only send if changed
        if (settingsStr !== this.state.lastSentSettings) {
            this.state.lastSentSettings = settingsStr;
            this.state.isProcessingSettings = true;

            Streamlit.setComponentValue({
                [StreamlitValue.SETTINGS]: settings,
                [StreamlitValue.SELECTION]: this.state.lastSentSelection,
            });

            setTimeout(() => {
                this.state.isProcessingSettings = false;
            }, 100);
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
        }

        // Structure onselect
        if (visualizer.structure) {
            this.state.originalStructOnselect = visualizer.structure.onselect;
            visualizer.structure.onselect = (indexes: any) => {
                this.state.originalStructOnselect?.(indexes);
                this.sendSelectionToStreamlit(indexes);
            };

            if (typeof visualizer.structure.activeChanged === 'function') {
                const originalActiveChanged = visualizer.structure.activeChanged.bind(
                    visualizer.structure
                );
                this.state.originalStructActiveChanged = originalActiveChanged;

                visualizer.structure.activeChanged = (guid: any, indexes: any) => {
                    originalActiveChanged?.(guid, indexes);
                    this.sendSelectionToStreamlit(indexes);
                };
            }
        }

        // Info onchange
        if (visualizer.info) {
            this.state.originalInfoOnchange = visualizer.info.onchange;
            visualizer.info.onchange = (indexes: any) => {
                this.state.originalInfoOnchange?.(indexes);
                this.sendSelectionToStreamlit(indexes);
            };
        }

        // Settings change - from visualizer to Streamlit
        visualizer.onSettingChange((_keys: string[], _value: unknown) => {
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
                if (selectedIndex !== null) {
                    // Small delay to ensure visualizer is fully loaded
                    setTimeout(() => {
                        this.applyExternalSelection(selectedIndex);
                    }, 500);
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
