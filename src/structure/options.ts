/**
 * @packageDocumentation
 * @module settings
 */

import assert from 'assert';

import { HTMLOption, OptionsGroup, SavedSettings } from '../options';
import { optionValidator } from '../options';
import { PositioningCallback, makeDraggable, sendWarning } from '../utils';

import BARS_SVG from '../static/bars.svg';
import HTML_OPTIONS from './options.html';

export class StructureOptions extends OptionsGroup {
    /// should we show bonds
    public bonds: HTMLOption<'boolean'>;
    /// should we use space filling representation
    public spaceFilling: HTMLOption<'boolean'>;
    /// should we show atoms labels
    public atomLabels: HTMLOption<'boolean'>;
    /// should we show unit cell information and lines
    public unitCell: HTMLOption<'boolean'>;
    /// number of repetitions in the `a/b/c` direction for the supercell
    public supercell: [HTMLOption<'int'>, HTMLOption<'int'>, HTMLOption<'int'>];
    // should we spin the representation
    public rotation: HTMLOption<'boolean'>;
    // which axis system to use (none, xyz, abc)
    public axes: HTMLOption<'string'>;
    // keep the orientation constant when loading a new structure
    public keepOrientation: HTMLOption<'boolean'>;
    // options related to environments
    public environments: {
        // should we display environments & environments options
        activated: HTMLOption<'boolean'>;
        // automatically center the environment when loading it
        center: HTMLOption<'boolean'>;
        // the cutoff value for spherical environments
        cutoff: HTMLOption<'number'>;
        // which style for atoms not in the environment
        bgStyle: HTMLOption<'string'>;
        // which colors for atoms not in the environment
        bgColor: HTMLOption<'string'>;
    };

    /// The HTML element containing the settings modal
    private _modal: HTMLElement;
    // Callback to get the initial positioning of the settings modal.
    private _positionSettingsModal: PositioningCallback;

    constructor(root: HTMLElement, guid: string, positionSettings: PositioningCallback) {
        super();

        this.bonds = new HTMLOption('boolean', true);
        this.spaceFilling = new HTMLOption('boolean', false);
        this.atomLabels = new HTMLOption('boolean', false);
        this.unitCell = new HTMLOption('boolean', false);
        this.rotation = new HTMLOption('boolean', false);

        this.supercell = [
            new HTMLOption('int', 1),
            new HTMLOption('int', 1),
            new HTMLOption('int', 1),
        ];

        const validateSupercell = (value: number) => {
            if (!(Number.isInteger(value) && value > 0)) {
                throw Error('supercell count should be a positive integer');
            }
        };
        this.supercell[0].validate = validateSupercell;
        this.supercell[1].validate = validateSupercell;
        this.supercell[2].validate = validateSupercell;

        this.axes = new HTMLOption('string', 'off');
        this.axes.validate = optionValidator(['off', 'abc', 'xyz'], 'axes');
        this.keepOrientation = new HTMLOption('boolean', false);

        this.environments = {
            activated: new HTMLOption('boolean', true),
            bgColor: new HTMLOption('string', 'grey'),
            bgStyle: new HTMLOption('string', 'ball-stick'),
            center: new HTMLOption('boolean', false),
            cutoff: new HTMLOption('number', 4.0),
        };

        this.environments.bgColor.validate = optionValidator(
            ['grey', 'CPK'],
            'background atoms coloring'
        );
        this.environments.bgStyle.validate = optionValidator(
            ['ball-stick', 'licorice', 'hide'],
            'background atoms style'
        );
        this.environments.cutoff.validate = (value) => {
            if (value < 0) {
                throw Error('cutoff should be a positive number');
            }
        };

        this._positionSettingsModal = positionSettings;

        this._modal = this._insertHTML(root, guid);
        document.body.appendChild(this._modal);
        this._bind(guid);
    }

    /**
     * Remove all HTML added by this [[StructureSettings]] in the current
     * document
     */
    public remove(): void {
        if (this._modal.classList.contains('show')) {
            const close = this._modal.querySelector('.close');
            assert(close !== null);
            (close as HTMLElement).click();
        }
        this._modal.remove();
    }

    /**
     * Applies saved settings, possibly filling in with default values
     */
    public applySettings(settings: SavedSettings): void {
        // don't warn for packedCell setting if is was set to false, which is
        // now the only possible way of doing it
        if ('packedCell' in settings) {
            if (settings.packedCell !== false) {
                sendWarning(
                    'packedCell option has been removed, but it is set to true in the settings'
                );
            }
            delete settings.packedCell;
        }

        super.applySettings(settings);
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
        const openSettings = template.content.firstChild as HTMLElement;
        root.append(openSettings);

        // replace id to ensure they are unique even if we have mulitple viewers
        // on a single page
        // prettier-ignore
        template.innerHTML = HTML_OPTIONS
            .replace(/id="(.*?)"/g, (_: string, id: string) => `id="${guid}-${id}"`)
            .replace(/for="(.*?)"/g, (_: string, id: string) => `for="${guid}-${id}"`)
            .replace(/data-target="#(.*?)"/g, (_: string, id: string) => `data-target="#${guid}-${id}"`);

        const modal = template.content.firstChild as HTMLElement;
        const modalDialog = modal.childNodes[1] as HTMLElement;
        assert(modalDialog !== undefined);
        assert(modalDialog.classList.contains('modal-dialog'));

        // Position modal near the actual viewer
        openSettings.addEventListener('click', () => {
            // only set style once, on first open, and keep previous position
            // on next open to keep the 'draged-to' position
            if (modalDialog.getAttribute('data-initial-modal-positions-set') === null) {
                modalDialog.setAttribute('data-initial-modal-positions-set', 'true');

                // display: block to ensure modalDialog.offsetWidth is non-zero
                (modalDialog.parentNode as HTMLElement).style.display = 'block';

                const { top, left } = this._positionSettingsModal(
                    modalDialog.getBoundingClientRect()
                );

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

    /** Bind all options to the corresponding HTML elements */
    private _bind(guid: string): void {
        this.atomLabels.bind(`${guid}-atom-labels`, 'checked');
        this.spaceFilling.bind(`${guid}-space-filling`, 'checked');
        this.bonds.bind(`${guid}-bonds`, 'checked');

        this.rotation.bind(`${guid}-rotation`, 'checked');
        this.unitCell.bind(`${guid}-unit-cell`, 'checked');

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
