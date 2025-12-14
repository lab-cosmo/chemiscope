import {
    Dataset,
    DefaultVisualizer,
    EnvironmentIndexer,
    GUID,
    Indexes,
    MapVisualizer,
    Settings,
    StructureVisualizer,
    Warnings,
} from '../../../src/index';

export type ChemiscopeMode = 'default' | 'structure' | 'map';

export interface ChemiscopeArgs {
    dataset: Dataset;
    height?: number;
    width?: number | string;
    selected_index?: number | undefined;
    mode?: ChemiscopeMode;
    settings?: Partial<Settings>;
    no_info_panel?: boolean;
}

export type SelectCallback = ((indexes: Indexes) => void) | null;
export type ActiveChangedCallback = ((guid: GUID, indexes: Indexes) => void) | null;

export interface ChemiscopeVisualizer {
    map?: {
        select: (indexes: Indexes) => void;
        onselect: SelectCallback;
        activeChanged?: ActiveChangedCallback;
    };
    structure?: {
        onselect: SelectCallback;
        activeChanged?: ActiveChangedCallback;
        activeIndex: number;
    };
    info?: {
        onchange: ((indexes: Indexes) => void) | null;
    };
    select: (indexes: Indexes) => void;
    applySettings: (settings: Partial<Settings>) => void;
    saveSettings: () => Partial<Settings>;
    onSettingChange: (callback: (keys: string[], value: unknown) => void) => void;
}

export interface VisualizerClass {
    load(
        config: unknown,
        dataset: Dataset,
        warnings: Warnings
    ): Promise<ChemiscopeVisualizer>;
}

export interface ChemiscopeGlobal {
    DefaultVisualizer: typeof DefaultVisualizer;
    StructureVisualizer: typeof StructureVisualizer;
    MapVisualizer: typeof MapVisualizer;
    EnvironmentIndexer: typeof EnvironmentIndexer;
    Warnings: typeof Warnings;
}

declare global {
    interface Window {
        Chemiscope?: ChemiscopeGlobal;
    }
}
