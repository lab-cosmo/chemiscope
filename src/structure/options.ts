/**
 * @packageDocumentation
 * @module settings
 */

import assert from 'assert';

import Collapse from '../collapse';
import Modal from '../modal';
import { Settings } from '../dataset';
import { HTMLOption, OptionsGroup } from '../options';
import { optionValidator } from '../options';
import { PositioningCallback, getByID, makeDraggable, sendWarning } from '../utils';

import BARS_SVG from '../static/bars.svg';
import HTML_OPTIONS from './options.html.in';

export class StructureOptions extends OptionsGroup {
    /// should we show bonds
    public bonds: HTMLOption<'boolean'>;
    /// should we show atoms
    public atoms: HTMLOption<'boolean'>;
    /// should we use space filling representation
    public spaceFilling: HTMLOption<'boolean'>;
    /// should we show atoms labels
    public atomLabels: HTMLOption<'boolean'>;
    /// if shapes are in the dataset, which set of shapes should we show
    public shape: HTMLOption<'string'>;
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
    // trajectory playback delay in seconds
    public playbackDelay: HTMLOption<'number'>;
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

    /// The Modal instance
    private _modal: Modal;
    /// The HTML element containing the button to open the settings modal
    private _openModal: HTMLElement;
    // Callback to get the initial positioning of the settings modal.
    private _positionSettingsModal: PositioningCallback;

    constructor(root: HTMLElement, positionSettings: PositioningCallback) {
        super();

        this.bonds = new HTMLOption('boolean', true);
        this.atoms = new HTMLOption('boolean', true);
        this.spaceFilling = new HTMLOption('boolean', false);
        this.atomLabels = new HTMLOption('boolean', false);
        this.shape = new HTMLOption('string', '');
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
        this.playbackDelay = new HTMLOption('number', 700);

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

        const { openModal, modal } = this._createSettingsHTML();
        this._modal = modal;
        this._modal.shadow.adoptedStyleSheets = (
            root.getRootNode() as ShadowRoot
        ).adoptedStyleSheets;
        this._openModal = openModal;
        root.appendChild(this._openModal);

        this._bind();
    }

    /** Get in a element in the modal from its id */
    public getModalElement<T extends HTMLElement = HTMLElement>(id: string): T {
        return getByID(id, this._modal.shadow);
    }

    /** Get the modal used for these settings */
    public get modal(): Modal {
        return this._modal;
    }

    /**
     * Remove all HTML added by this {@link StructureSettings} in the current
     * document
     */
    public remove(): void {
        this._modal.remove();
        this._openModal.remove();
    }

    /**
     * Applies saved settings, possibly filling in with default values
     */
    public applySettings(settings: Settings): void {
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
     * Create HTML needed for structure settings.
     *
     * The HTML elements are returned, not yet inserted in the document.
     *
     * @return the HTML element containing the setting modal, and the button to open the modal
     */
    private _createSettingsHTML(): { modal: Modal; openModal: HTMLElement } {
        // use HTML5 template to generate a DOM object from an HTML string
        const template = document.createElement('template');
        template.innerHTML = `<button
            class="btn btn-light btn-sm chsp-viewer-button"
            style="top: 4px; right: 4px; opacity: 1;">
                <div>${BARS_SVG}</div>
            </button>`;
        const openModal = template.content.firstChild as HTMLElement;

        template.innerHTML = HTML_OPTIONS;
        const modalElement = template.content.querySelector('.modal');
        assert(modalElement !== null && modalElement instanceof HTMLElement);
        const modalDialog = modalElement.querySelector('.modal-dialog');
        assert(modalDialog !== null && modalDialog instanceof HTMLElement);

        // Position modal near the actual viewer
        openModal.addEventListener('click', () => {
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

            modal.open();
        });

        // make the settings modal draggable
        makeDraggable(modalDialog, '.modal-header');

        // Stop propagation of keydown events. This is required for the Jupyter integration,
        // otherwise jupyter tries to interpret key press in the modal as its own input
        modalElement.addEventListener('keydown', (event) => event.stopPropagation());

        const modal = new Modal(modalElement);
        Collapse.initialize(modalElement);

        return { modal, openModal };
    }

    /** Bind all options to the corresponding HTML elements */
    private _bind(): void {
        this.atomLabels.bind(this.getModalElement('atom-labels'), 'checked');

        const selectShape = this.getModalElement<HTMLSelectElement>('shapes');
        this.shape.bind(selectShape, 'options');

        this.spaceFilling.bind(this.getModalElement('space-filling'), 'checked');
        this.bonds.bind(this.getModalElement('bonds'), 'checked');
        this.atoms.bind(this.getModalElement('atoms'), 'checked');

        this.rotation.bind(this.getModalElement('rotation'), 'checked');
        this.unitCell.bind(this.getModalElement('unit-cell'), 'checked');

        this.supercell[0].bind(this.getModalElement('supercell-a'), 'value');
        this.supercell[1].bind(this.getModalElement('supercell-b'), 'value');
        this.supercell[2].bind(this.getModalElement('supercell-c'), 'value');

        this.axes.bind(this.getModalElement('axes'), 'value');
        this.keepOrientation.bind(this.getModalElement('keep-orientation'), 'checked');
        this.playbackDelay.bind(this.getModalElement('playback-delay'), 'value');

        this.environments.activated.bind(this.getModalElement('env-activated'), 'checked');
        this.environments.bgColor.bind(this.getModalElement('env-bg-color'), 'value');
        this.environments.bgStyle.bind(this.getModalElement('env-bg-style'), 'value');
        this.environments.cutoff.bind(this.getModalElement('env-cutoff'), 'value');
        this.environments.center.bind(this.getModalElement('env-center'), 'checked');
    }
}
