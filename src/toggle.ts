import { getElement } from './utils';
import * as styles from './styles';

export class DisplayToggle {
    /// Shadow root for isolation
    private _shadow: ShadowRoot;
    /// Toggle input element
    private _toggle: HTMLInputElement;
    /// Callback fired when the user changes the toggle value
    public onchange: (checked: boolean) => void;

    /**
     * Create a new {@link DisplayToggle} instance
     *
     * @param element HTML element or HTML id of the DOM element where the toggle will be attached
     * @param toggled flag indicating if the toggle should be checked
     * @param disabled flag indicating if the toggle should be disabled
     */
    constructor(
        element: string | HTMLElement,
        toggled: boolean = false,
        disabled: boolean = false
    ) {
        // Create a container element
        const containerElement = getElement(element);
        const hostElement = document.createElement('div');
        containerElement.appendChild(hostElement);

        // Attach to shadow
        this._shadow = hostElement.attachShadow({ mode: 'open' });
        this._shadow.adoptedStyleSheets = [styles.bootstrap, styles.chemiscope];

        // Create a toggle element
        this._toggle = this._createToggleElement(toggled, disabled);
        this._shadow.appendChild(this._toggle.parentElement as Node);

        // Set up events
        this._toggle.onchange = () => this.onchange?.(this._toggle.checked);
        this.onchange = () => {};
    }

    /**
     * Create the HTML structure for the toggle element
     *
     * @param toggled flag indicating if the toggle should be checked
     * @param disabled flag indicating if the toggle should be disabled
     * @returns the input element of the toggle
     */
    private _createToggleElement(toggled: boolean, disabled: boolean): HTMLInputElement {
        const toggleElement = document.createElement('div');
        toggleElement.innerHTML = `
            <div class="chsp-mode-toggle form-check form-switch">
                <input class="form-check-input" id="atom-display-toggle" type="checkbox" 
                       ${disabled ? 'disabled' : ''} ${toggled ? 'checked' : ''} />
                <label class="form-check-label" for="atom-display-toggle" title="atom display">
                    atom display
                </label>
            </div>
        `;
        return toggleElement.querySelector('input') as HTMLInputElement;
    }

    /**
     * Remove HTML added by DisplayToggle in the current document
     */
    remove(): void {
        this._shadow.host.remove();
    }
}
