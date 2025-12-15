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
import { Dataset, EnvironmentIndexer, GUID, Indexes, Settings } from 'chemiscope';
import {
    ActiveChangedCallback,
    ChemiscopeArgs,
    ChemiscopeGlobal,
    ChemiscopeMode,
    ChemiscopeVisualizer,
    SelectCallback,
    VisualizerClass,
} from './types';

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
        if (selectedIndex === undefined) {
            return;
        }

        // Only process if selection actually changed
        if (selectedIndex !== this.state.currentSelection) {
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
        if (!settings || !this.state.visualizer) {
            return;
        }

        const settingsStr = JSON.stringify(settings);

        // Only process if settings actually changed
        if (settingsStr !== this.state.currentSettings) {
            this.state.settingsFromPython = true;
            try {
                this.state.currentSettings = settingsStr;
                const settingsToApply = this.convertStructureSettingsToArray(settings);
                this.state.visualizer.applySettings(settingsToApply);
            } finally {
                this.state.settingsFromPython = false;
            }
        }
    }

    private convertStructureSettingsToArray(settings: Partial<Settings>): Partial<Settings> {
        // if structure is already an array or not present
        if (!settings.structure || Array.isArray(settings.structure)) {
            return settings;
        }

        const currentSettings = this.state.visualizer?.saveSettings() || {};
        const currentStructures: Settings[] = Array.isArray(currentSettings.structure)
            ? (currentSettings.structure as Settings[])
            : [];

        if (currentStructures.length === 0) {
            return {
                ...settings,
                structure: [settings.structure as Settings],
            };
        }

        const activeIndex = this.getActiveStructureIndex();
        currentStructures[activeIndex] = settings.structure as Settings;

        return {
            ...settings,
            structure: currentStructures,
        };
    }

    private applySelection(selectedIndex: number | null | undefined): void {
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

        this.state.originalSelect?.(indexes);
    }

    private sendSelectionToStreamlit(indexes: Indexes | null): void {
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

            // Get current settings with only the active structure's settings
            const settingsToSend = this.getSettingsWithActiveStructure();
            const settingsStr = JSON.stringify(settingsToSend);
            this.state.currentSettings = settingsStr;

            Streamlit.setComponentValue({
                [StreamlitValue.SETTINGS]: settingsToSend,
                [StreamlitValue.SELECTION]: structureIdToSend,
            });
        }
    }

    private sendSettingsToStreamlit(): void {
        if (this.state.settingsFromPython) {
            return;
        }

        const settingsToSend = this.getSettingsWithActiveStructure();
        const settingsStr = JSON.stringify(settingsToSend);

        // Only send if changed
        if (settingsStr !== this.state.currentSettings) {
            this.state.currentSettings = settingsStr;

            Streamlit.setComponentValue({
                [StreamlitValue.SETTINGS]: settingsToSend,
                [StreamlitValue.SELECTION]: this.state.currentSelection,
            });
        }
    }

    private getSettingsWithActiveStructure(): Settings {
        const currentSettings = this.state.visualizer?.saveSettings() || {};
        const activeIndex = this.getActiveStructureIndex();

        // if there are structure settings, extract only the active one
        if (Array.isArray(currentSettings.structure) && currentSettings.structure.length > 0) {
            const activeStructureSettings =
                currentSettings.structure[activeIndex] || currentSettings.structure[0];
            return {
                ...currentSettings,
                structure: activeStructureSettings,
            };
        }

        return currentSettings;
    }

    private getActiveStructureIndex(): number {
        const visualizer = this.state.visualizer;
        if (!visualizer?.structure) {
            return 0;
        }

        return visualizer.structure.activeIndex;
    }

    private setupBidirectionalSync(): void {
        const visualizer = this.state.visualizer;
        if (!visualizer) {
            return;
        }

        // Map onselect
        if (visualizer.map) {
            const originalMapOnselect = visualizer.map.onselect?.bind(visualizer.map) ?? null;

            this.state.originalMapOnselect = originalMapOnselect;
            visualizer.map.onselect = (indexes: Indexes) => {
                originalMapOnselect?.(indexes);
                this.sendSelectionToStreamlit(indexes);
            };

            if (typeof visualizer.map.activeChanged === 'function') {
                const originalActiveChanged = visualizer.map.activeChanged.bind(visualizer.map);
                this.state.originalMapActiveChanged = originalActiveChanged;

                visualizer.map.activeChanged = (guid: GUID, indexes: Indexes) => {
                    originalActiveChanged?.(guid, indexes);
                    this.sendSettingsToStreamlit();
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
                originalStructOnselect?.(indexes);
                this.sendSelectionToStreamlit(indexes);
            };

            if (typeof visualizer.structure.activeChanged === 'function') {
                const originalActiveChanged = visualizer.structure.activeChanged.bind(
                    visualizer.structure
                );
                this.state.originalStructActiveChanged = originalActiveChanged;

                visualizer.structure.activeChanged = (guid: GUID, indexes: Indexes) => {
                    originalActiveChanged?.(guid, indexes);
                    this.sendSettingsToStreamlit();
                    this.sendSelectionToStreamlit(indexes);
                };
            }
        }

        // Info onchange
        if (visualizer.info) {
            const originalSelect = visualizer.select.bind(visualizer);
            this.state.originalSelect = originalSelect;

            visualizer.select = (indexes: Indexes) => {
                originalSelect(indexes);
                this.sendSelectionToStreamlit(indexes);
            };
        }

        // Settings change - from visualizer to Streamlit
        visualizer.onSettingChange(() => {
            this.sendSettingsToStreamlit();
        });
    }

    private initializeVisualizer(
        Chemiscope: ChemiscopeGlobal,
        dataset: Dataset,
        mode: ChemiscopeMode,
        settings: Partial<Settings>,
        selectedIndex: number | null | undefined
    ): void {
        // Merge initial settings
        try {
            dataset.settings = Object.assign({}, dataset.settings, settings);
        } catch (e) {
            // eslint-disable-next-line no-console
            console.warn('Could not attach settings to dataset:', String(e));
        }

        const warnings = new Chemiscope.Warnings();
        warnings.addHandler((message: string, timeout: number = 4000) => {
            displayWarning(message, timeout);
        });

        toggleLoadingVisible(true);

        const visualizerClass = getVisualizerClass(mode, Chemiscope);

        // eslint-disable-next-line @typescript-eslint/no-unsafe-call
        visualizerClass
            .load(CONFIG, dataset, warnings)
            .then((v: ChemiscopeVisualizer) => {
                this.state.visualizer = v;
                this.state.indexer = new Chemiscope.EnvironmentIndexer(
                    dataset.structures,
                    dataset.environments
                );

                this.setupBidirectionalSync();

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
                const message = err instanceof Error ? err.message : String(err);

                // eslint-disable-next-line no-console
                console.error('Error loading visualizer:', message);
                displayWarning(`Error loading visualization: ${message}`);
            })
            .finally(() => {
                toggleLoadingVisible(false);
            });
    }
}
