/**
 * @packageDocumentation
 * @module settings
 */

import {HTMLSetting, SettingsGroup, SettingsPreset, settingsValidator} from '../settings';
import {makeDraggable, PositioningCallback} from '../utils';

import BARS_SVG from '../static/bars.svg';
import HTML_SETTINGS from './settings.html';

export class StructureSettings extends SettingsGroup {
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
        presets: SettingsPreset = {},
    ) {
        super();

        this.bonds = new HTMLSetting('boolean', true);
        this.spaceFilling = new HTMLSetting('boolean', false);
        this.atomLabels = new HTMLSetting('boolean', false);
        this.unitCell = new HTMLSetting('boolean', false);
        this.packedCell = new HTMLSetting('boolean', false);
        this.rotation = new HTMLSetting('boolean', false);

        this.supercell = [
            new HTMLSetting('int', 1),
            new HTMLSetting('int', 1),
            new HTMLSetting('int', 1),
        ];

        const validateSupercell = (value: number) => {
            if (!(Number.isInteger(value) && value > 0)) {
                throw Error('supercell count should be a positive integer');
            }
        };
        this.supercell[0].validate = validateSupercell;
        this.supercell[1].validate = validateSupercell;
        this.supercell[2].validate = validateSupercell;

        this.axes = new HTMLSetting('string', 'off');
        this.axes.validate = settingsValidator(['off', 'abc', 'xyz'], 'axes');
        this.keepOrientation = new HTMLSetting('boolean', false);

        this.environments = {
            activated: new HTMLSetting('boolean', true),
            bgColor:  new HTMLSetting('string', 'grey'),
            bgStyle: new HTMLSetting('string', 'licorice'),
            center: new HTMLSetting('boolean', false),
            cutoff: new HTMLSetting('number', 4.0),
        };

        this.environments.bgColor.validate = settingsValidator(
            ['grey', 'CPK'], 'background atoms coloring',
        );
        this.environments.bgStyle.validate = settingsValidator(
            ['licorice', 'ball-stick', 'hide'], 'background atoms style',
        );
        this.environments.cutoff.validate = (value) => {
            if (value < 0) {
                throw Error('cutoff should be a positive number');
            }
        };

        this._positionSettingsModal = positionSettings;

        this._settingsModal = this._insertHTML(root, guid);
        document.body.appendChild(this._settingsModal);
        this._bindSettings(guid);

        this.applyPresets(presets);
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
