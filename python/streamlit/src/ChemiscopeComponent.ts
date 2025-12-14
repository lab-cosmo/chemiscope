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
import { Dataset, EnvironmentIndexer, GUID, Indexes, Settings } from '../../../src/index';
import {
    ActiveChangedCallback,
    ChemiscopeArgs,
    ChemiscopeGlobal,
    ChemiscopeMode,
    ChemiscopeVisualizer,
    SelectCallback,
    VisualizerClass,
} from '../types/chemiscope';

enum StreamlitValue {
    SELECTION = 'selected_id',
    SETTINGS = 'settings',
}

function getChemiscope(): ChemiscopeGlobal | null {
    if (!window.Chemiscope) {
        // eslint-disable-next-line no-console
        console.error('window.Chemiscope not found.');
        return null;
    }
    return window.Chemiscope;
}

function getVisualizerClass(mode: ChemiscopeMode, Chemiscope: ChemiscopeGlobal): VisualizerClass {
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
        indexer: null as EnvironmentIndexer | null,
        loaded: false,

        // Track selections
        currentSelection: null as number | null,
        currentActive: null as string | null,

        // Track settings
        currentSettings: null as string | null,

        // Track source of the selection events to avoid loops
        selectionFromPython: false,
        settingsFromPython: false,

        // Original callbacks
        originalMapOnselect: null as SelectCallback,
        originalMapActiveChanged: null as ActiveChangedCallback,
        originalStructOnselect: null as SelectCallback,
        originalStructActiveChanged: null as ActiveChangedCallback,
        originalSelect: null as SelectCallback,
    };

    constructor() {
        Streamlit.events.addEventListener(Streamlit.RENDER_EVENT, (event: Event) => {
            this.onRender((event as CustomEvent<RenderData>).detail);
        });
        Streamlit.setComponentReady();
    }

    private onRender(data: RenderData): void {
        // eslint-disable-next-line no-console
        console.log('********* RENDERING STREAMLIT COMPONENT ************');
        const args = data.args as ChemiscopeArgs;
        const dataset = args.dataset;

        if (!dataset) {
            return;
        }

        const mode = typeof args.mode === 'string' ? args.mode : 'default';
        const noInfo = typeof args.no_info_panel === 'boolean' ? args.no_info_panel : false;
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
                heightArg,
                noInfo
            );
        } else {
            this.handleUpdate(args, widthArg, heightArg);
        }
    }

    private handleFirstRender(
        Chemiscope: ChemiscopeGlobal,
        dataset: Dataset,
        mode: ChemiscopeMode,
        settings: Partial<Settings>,
        selectedIndex: number | undefined,
        widthArg: string | number,
        heightArg: number,
        noInfo: boolean
    ): void {
        this.state.loaded = true;

        const root = getOrCreateRoot();
        root.innerHTML = generateHTMLForMode(mode, noInfo);

        applyWidthPolicy(widthArg, root);
        applyHeightPolicy(heightArg, root);

        // Store initial received values
        this.state.currentSelection = selectedIndex !== undefined ? selectedIndex : null;
        this.state.currentSettings = JSON.stringify(settings);

        this.initializeVisualizer(Chemiscope, dataset, mode, settings, selectedIndex);
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
        // eslint-disable-next-line no-console
        console.log('Handling selection update:', selectedIndex);
        if (selectedIndex === undefined) {
            return;
        }

        // Only process if selection actually changed
        if (selectedIndex !== this.state.currentSelection) {
            // eslint-disable-next-line no-console
            console.log('update selection');

            this.state.selectionFromPython = true;
            try {
                this.state.currentSelection = selectedIndex;
                this.applySelection(selectedIndex);
            } finally {
                this.state.selectionFromPython = false;
            }
        }
    }

    private handleSettingsUpdate(settings: Partial<Settings> | undefined): void {
        // eslint-disable-next-line no-console
        console.log('Handling settings update');
        if (!settings || !this.state.visualizer) {
            return;
        }

        const settingsStr = JSON.stringify(settings);

        // Only process if settings actually changed
        if (settingsStr !== this.state.currentSettings) {
            // eslint-disable-next-line no-console
            console.log('update settings');

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
        // eslint-disable-next-line no-console
        console.log(
            'Applying selection:',
            selectedIndex,
            'from python',
            this.state.selectionFromPython
        );
        const { visualizer, indexer } = this.state;
        if (!visualizer || !indexer) {
            return;
        }

        if (selectedIndex === null || selectedIndex === undefined) {
            // Handle deselection if needed
            return;
        }

        const indexes = indexer.fromStructureAtom('structure', selectedIndex);
        if (!indexes) {
            // eslint-disable-next-line no-console
            console.warn('No environment for structure index', selectedIndex);
            return;
        }

        // eslint-disable-next-line no-console
        console.log('calling originalselect, ', indexes);
        this.state.originalSelect?.(indexes);
    }

    private sendSelectionToStreamlit(indexes: Indexes | null): void {
        // eslint-disable-next-line no-console
        console.log('Sending selection to Streamlit:', indexes, this.state.selectionFromPython);
        if (this.state.selectionFromPython) {
            return;
        }

        let structureIdToSend: number | null = null;
        if (indexes && typeof indexes.structure === 'number') {
            structureIdToSend = indexes.structure;
        }

        // Only send if changed
        // eslint-disable-next-line no-console
        console.log('send:', structureIdToSend, 'current ', this.state.currentSelection);
        if (structureIdToSend !== this.state.currentSelection) {
            this.state.currentSelection = structureIdToSend;

            // Get current settings
            const currentSettings = this.state.visualizer?.saveSettings() || {};
            const settingsStr = JSON.stringify(currentSettings);
            this.state.currentSettings = settingsStr;

            Streamlit.setComponentValue({
                [StreamlitValue.SETTINGS]: currentSettings,
                [StreamlitValue.SELECTION]: structureIdToSend,
            });
        }
    }

    private sendSettingsToStreamlit(settings: Settings): void {
        // eslint-disable-next-line no-console
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
            const originalMapOnselect = visualizer.map.onselect?.bind(visualizer.map) ?? null;

            this.state.originalMapOnselect = originalMapOnselect;
            visualizer.map.onselect = (indexes: Indexes) => {
                // eslint-disable-next-line no-console
                console.log('map.onselect');
                originalMapOnselect?.(indexes);
                this.sendSelectionToStreamlit(indexes);
            };

            if (typeof visualizer.map.activeChanged === 'function') {
                const originalActiveChanged = visualizer.map.activeChanged.bind(visualizer.map);
                this.state.originalMapActiveChanged = originalActiveChanged;

                visualizer.map.activeChanged = (guid: GUID, indexes: Indexes) => {
                    // eslint-disable-next-line no-console
                    console.log('map.activeChanged', guid, indexes);
                    originalActiveChanged?.(guid, indexes);
                    this.sendSelectionToStreamlit(indexes);
                };
            }
        }

        // Structure onselect - selection from structure viewer to Streamlit
        if (visualizer.structure) {
            const originalStructOnselect =
                visualizer.structure.onselect?.bind(visualizer.structure) ?? null;
            this.state.originalStructOnselect = originalStructOnselect;
            visualizer.structure.onselect = (indexes: Indexes) => {
                // eslint-disable-next-line no-console
                console.log('structure.onselect');
                originalStructOnselect?.(indexes);
                this.sendSelectionToStreamlit(indexes);
            };

            if (typeof visualizer.structure.activeChanged === 'function') {
                const originalActiveChanged = visualizer.structure.activeChanged.bind(
                    visualizer.structure
                );
                this.state.originalStructActiveChanged = originalActiveChanged;

                visualizer.structure.activeChanged = (guid: GUID, indexes: Indexes) => {
                    // eslint-disable-next-line no-console
                    console.log('structure.activeChanged', guid, indexes);
                    originalActiveChanged?.(guid, indexes);
                    this.sendSelectionToStreamlit(indexes);
                };
            }
        }

        // Info onchange
        if (visualizer.info) {
            const originalSelect = visualizer.select.bind(visualizer);
            this.state.originalSelect = originalSelect;

            visualizer.select = (indexes: Indexes) => {
                // eslint-disable-next-line no-console
                console.log('visualizer.select');

                originalSelect(indexes);
                this.sendSelectionToStreamlit(indexes);
            };
        }

        // Settings change - from visualizer to Streamlit
        visualizer.onSettingChange(() => {
            // eslint-disable-next-line no-console
            console.log('onsettingschange');
            const currentSettings = visualizer.saveSettings();
            this.sendSettingsToStreamlit(currentSettings);
        });
    }

    private initializeVisualizer(
        Chemiscope: ChemiscopeGlobal,
        dataset: Dataset,
        mode: ChemiscopeMode,
        settings: Partial<Settings>,
        selectedIndex: number | null | undefined
    ): void {
        const visualizerClass = getVisualizerClass(mode, Chemiscope);

        // Merge initial settings
        try {
            dataset.settings = Object.assign({}, dataset.settings, settings);
        } catch (e) {
            // eslint-disable-next-line no-console
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
                // eslint-disable-next-line no-console
                console.error('Error loading visualizer:', err);
                displayWarning('Error loading visualization: ' + String(err));
            })
            .finally(() => {
                toggleLoadingVisible(false);
            });
    }
}
