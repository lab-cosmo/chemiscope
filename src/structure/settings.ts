/**
 * @packageDocumentation
 * @module settings
 */

import {HTMLSetting, PositioningCallback} from '../utils';
import {makeDraggable} from '../utils';

import BARS_SVG from '../static/bars.svg';
import HTML_SETTINGS from './settings.html';

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

    /// The HTML element containing the settings modal
    private _settingsModal: HTMLElement;
    // Callback to get the initial positioning of the settings modal.
    private _positionSettingsModal: PositioningCallback;

    constructor(
        root: HTMLElement,
        guid: string,
        positionSettings: PositioningCallback,
        presets: Partial<StructurePresets> = {},
    ) {
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

        this._positionSettingsModal = positionSettings;

        this._settingsModal = this._insertHTML(root, guid);
        document.body.appendChild(this._settingsModal);
        this._bindSettings(guid);

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

    /**
     * Remove all HTML added by this [[StructureSettings]] in the current
     * document
     */
    public remove(): void {
        this._settingsModal.remove();
    }

    /**
     * Insert HTML needed for setting in the page, adding the "open settings"
     * button to `root`. The setting modal HTML element is returned, not yet
     * inserted in the document.
     *
     * @param  root where to place the HTML button
     * @param  guid unique identifier of the corresponding JSmolWidget, used as
     *              prefix for all elements ID
     * @return      the HTML element containing the setting modal
     */
    private _insertHTML(root: HTMLElement, guid: string): HTMLElement {
        // use HTML5 template to generate a DOM object from an HTML string
        const template = document.createElement('template');
        template.innerHTML = `<button
            class="btn btn-light btn-sm chsp-viewer-button"
            data-target="#${guid}-settings"
            data-toggle="modal"
            style="top: 5px; right: 5px; opacity: 1;">
                <div>${BARS_SVG}</div>
            </button>`;
        const openSettings = template.content.firstChild!;
        root.append(openSettings);

        // replace id to ensure they are unique even if we have mulitple viewers
        // on a single page
        template.innerHTML = HTML_SETTINGS
            .replace(/id=(.*?) /g, (_: string, id: string) => `id=${guid}-${id} `)
            .replace(/for=(.*?) /g, (_: string, id: string) => `for=${guid}-${id} `)
            .replace(/data-target=#(.*?) /g, (_: string, id: string) => `data-target=#${guid}-${id} `);

        const modal = template.content.firstChild! as HTMLElement;

        const modalDialog = modal.childNodes[1]! as HTMLElement;
        if (!modalDialog.classList.contains('modal-dialog')) {
            throw Error('internal error: missing modal-dialog class');
        }

        // Position modal near the actual viewer
        openSettings.addEventListener('click', () => {
            // only set style once, on first open, and keep previous position
            // on next open to keep the 'draged-to' position
            if (modalDialog.getAttribute('data-initial-modal-positions-set') === null) {
                modalDialog.setAttribute('data-initial-modal-positions-set', 'true');

                // display: block to ensure modalDialog.offsetWidth is non-zero
                (modalDialog.parentNode as HTMLElement).style.display = 'block';

                const {top, left} = this._positionSettingsModal(modalDialog.getBoundingClientRect());

                // set width first, since setting position can influence it
                modalDialog.style.width = `${modalDialog.offsetWidth}px`;
                // unset margins when using position: fixed
                modalDialog.style.margin = '0';
                modalDialog.style.position = 'fixed';
                modalDialog.style.top = `${top}px`;
                modalDialog.style.left = `${left}px`;
            }
        });

        // make the settings modal draggable
        makeDraggable(modalDialog, '.modal-header');

        return modal;
    }

    private _bindSettings(guid: string): void {
        // bind all the settings to corresponding HTML elements
        this.atomLabels.bind(`${guid}-atom-labels`, 'checked');
        this.spaceFilling.bind(`${guid}-space-filling`, 'checked');
        this.bonds.bind(`${guid}-bonds`, 'checked');

        this.rotation.bind(`${guid}-rotation`, 'checked');
        this.unitCell.bind(`${guid}-unit-cell`, 'checked');
        this.packedCell.bind(`${guid}-packed-cell`, 'checked');

        this.supercell[0].bind(`${guid}-supercell-a`, 'value');
        this.supercell[1].bind(`${guid}-supercell-b`, 'value');
        this.supercell[2].bind(`${guid}-supercell-c`, 'value');

        this.axes.bind(`${guid}-axes`, 'value');
        this.keepOrientation.bind(`${guid}-keep-orientation`, 'checked');

        this.environments.activated.bind(`${guid}-env-activated`, 'checked');
        this.environments.bgColor.bind(`${guid}-env-bg-color`, 'value');
        this.environments.bgStyle.bind(`${guid}-env-bg-style`, 'value');
        this.environments.cutoff.bind(`${guid}-env-cutoff`, 'value');
        this.environments.center.bind(`${guid}-env-center`, 'checked');
    }
}
