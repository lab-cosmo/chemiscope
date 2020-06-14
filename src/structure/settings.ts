/**
 * @packageDocumentation
 * @module settings
 */

import {HTMLSetting} from '../utils';

interface EnvironmentPresets {
    activated: boolean;
    center: boolean;
    cutoff: number;
    bgStyle: string;
    bgColor: string;
}

export interface StructurePresets {
    bonds: boolean;
    spaceFilling: boolean;
    atomLabels: boolean;
    unitCell: boolean;
    packedCell: boolean;
    supercell: number[];
    rotation: boolean;
    axes: string;
    environments: Partial<EnvironmentPresets>;
    keepOrientation: boolean;
}

export const STRUCTURE_DEFAULTS: StructurePresets = {
    atomLabels: false,
    axes: 'off',
    bonds: true,
    environments: {
        activated: true,
        bgColor: 'grey',
        bgStyle: 'licorice',
        center: false,
        cutoff: 4.0,
    },
    keepOrientation: false,
    packedCell: false,
    rotation: false,
    spaceFilling: false,
    supercell: [1, 1, 1],
    unitCell: false,
};

export class StructureSettings {
    // should we show bonds
    public bonds: HTMLSetting<'boolean'>;
    // should we use space filling representation
    public spaceFilling: HTMLSetting<'boolean'>;
    // should we show atoms labels
    public atomLabels: HTMLSetting<'boolean'>;
    // should we show unit cell information and lines
    public unitCell: HTMLSetting<'boolean'>;
    /// Is the current unit cell displayed as a packed cell?
    public packedCell: HTMLSetting<'boolean'>;
    /// number of repetitions in the `a/b/c` direction for the supercell
    public supercell: [HTMLSetting<'int'>, HTMLSetting<'int'>, HTMLSetting<'int'>];
    // should we spin the represenation
    public rotation: HTMLSetting<'boolean'>;
    // which axis system to use (none, xyz, abc)
    public axes: HTMLSetting<'string'>;
    // keep the orientation constant when loading a new structure if checked
    public keepOrientation: HTMLSetting<'boolean'>;
    // options related to environments
    public environments: {
        // should we display environments & environments options
        activated: HTMLSetting<'boolean'>;
        // automatically center the environment when loading it
        center: HTMLSetting<'boolean'>;
        // the cutoff value for spherical environments
        cutoff: HTMLSetting<'number'>;
        // which style for atoms not in the environment
        bgStyle: HTMLSetting<'string'>;
        // which colors for atoms not in the environment
        bgColor: HTMLSetting<'string'>;
    };

    constructor(presets: Partial<StructurePresets> = {}) {
        this.bonds = new HTMLSetting('boolean', STRUCTURE_DEFAULTS.bonds);
        this.spaceFilling = new HTMLSetting('boolean', STRUCTURE_DEFAULTS.spaceFilling);
        this.atomLabels = new HTMLSetting('boolean', STRUCTURE_DEFAULTS.atomLabels);
        this.unitCell = new HTMLSetting('boolean', STRUCTURE_DEFAULTS.unitCell);
        this.packedCell = new HTMLSetting('boolean', STRUCTURE_DEFAULTS.packedCell);
        this.supercell = [
            new HTMLSetting('int', STRUCTURE_DEFAULTS.supercell[0]),
            new HTMLSetting('int', STRUCTURE_DEFAULTS.supercell[1]),
            new HTMLSetting('int', STRUCTURE_DEFAULTS.supercell[2]),
        ];

        this.rotation = new HTMLSetting('boolean', STRUCTURE_DEFAULTS.rotation);
        this.axes = new HTMLSetting('string', STRUCTURE_DEFAULTS.axes);
        this.keepOrientation = new HTMLSetting('boolean', STRUCTURE_DEFAULTS.keepOrientation);

        const ENVIRONMENTS_DEFAUT = STRUCTURE_DEFAULTS.environments as EnvironmentPresets;
        this.environments = {
            activated: new HTMLSetting('boolean', ENVIRONMENTS_DEFAUT.activated),
            bgColor:  new HTMLSetting('string', ENVIRONMENTS_DEFAUT.bgColor),
            bgStyle: new HTMLSetting('string', ENVIRONMENTS_DEFAUT.bgStyle),
            center: new HTMLSetting('boolean', ENVIRONMENTS_DEFAUT.center),
            cutoff: new HTMLSetting('number', ENVIRONMENTS_DEFAUT.cutoff),
        };

        this.applyPresets(presets);
    }

    /**
     * Applies presets, possibly filling in with default values
     */
    public applyPresets(presets: Partial<StructurePresets> = {}) {
        const initial: StructurePresets = {
            ...STRUCTURE_DEFAULTS,
            ...presets,
        };
        // also complete the "environments" section
        initial.environments = {
            ...STRUCTURE_DEFAULTS.environments,
            ...initial.environments,
        };

        this.bonds.value = initial.bonds;
        this.spaceFilling.value = initial.spaceFilling;
        this.atomLabels.value = initial.atomLabels;
        this.unitCell.value = initial.unitCell;
        this.packedCell.value = initial.packedCell;
        this.supercell[0].value = initial.supercell[0];
        this.supercell[1].value = initial.supercell[1];
        this.supercell[2].value = initial.supercell[2];
        this.rotation.value = initial.rotation;
        this.axes.value = initial.axes;
        this.keepOrientation.value = initial.keepOrientation;
        this.environments.activated.value = initial.environments.activated!;
        this.environments.center.value = initial.environments.center!;
        this.environments.cutoff.value = initial.environments.cutoff!;
        this.environments.bgStyle.value = initial.environments.bgStyle!;
        this.environments.bgColor.value = initial.environments.bgColor!;
    }

    /**
     * Dumps presets, in a way that can e.g. be serialized to json
     */
    public dumpPresets(): StructurePresets {
        return {
            atomLabels: this.atomLabels.value,
            axes: this.axes.value,
            bonds: this.bonds.value,
            environments: {
                activated: this.environments.activated.value,
                bgColor: this.environments.bgColor.value,
                bgStyle: this.environments.bgStyle.value,
                center: this.environments.center.value,
                cutoff: this.environments.cutoff.value,
            },
            keepOrientation: this.keepOrientation.value,
            packedCell: this.packedCell.value,
            rotation: this.rotation.value,
            spaceFilling: this.spaceFilling.value,
            supercell: [ this.supercell[0].value, this.supercell[1].value, this.supercell[2].value],
            unitCell: this.unitCell.value,
        };
    }
}
