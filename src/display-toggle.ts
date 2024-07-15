import { getElement } from './utils';
import * as styles from './styles';
import { DisplayMode } from './indexer';

export class DisplayToggle {
    /// Shadow root for isolation
    private _shadow: ShadowRoot;
    /// Toggle buttons container element
    private _toggleContainer: HTMLElement;
    /// Callback fired when the user changes the toggle value
    public onchange: (mode: DisplayMode) => void;

    /**
     * Create a new {@link DisplayToggle} instance
     *
     * @param element HTML element or HTML id of the DOM element where the toggle will be attached
     * @param isPerAtom flag indicating if the atom mode should be checked
     */
    constructor(element: string | HTMLElement, isPerAtom: boolean = true) {
        // Create a container element
        const containerElement = getElement(element);
        const hostElement = document.createElement('div');
        containerElement.appendChild(hostElement);

        // Attach to shadow
        this._shadow = hostElement.attachShadow({ mode: 'open' });
        this._shadow.adoptedStyleSheets = [styles.bootstrap, styles.chemiscope];

        // Create a toggle element
        this._toggleContainer = this._createToggleElement(isPerAtom);
        this._shadow.appendChild(this._toggleContainer);

        // Set up events
        this.onchange = () => {};
    }

    /**
     * Create the HTML structure for the toggle element
     *
     * @param isPerAtom flag indicating if the atom mode should be checked
     * @returns the container element of the toggle buttons
     */
    private _createToggleElement(isPerAtom: boolean): HTMLElement {
        const toggleContainer = document.createElement('div');
        toggleContainer.innerHTML = `
            <div class="chsp-mode-toggle">
                <!-- Buttons -->
                <div class="btn-group-sm" role="group" aria-label="Mode toggle">
                    <button type="button" class="btn btn-outline-secondary ${!isPerAtom ? 'active' : ''}" id="structure-btn">Structures</button>
                    <button type="button" class="btn btn-outline-secondary ${isPerAtom ? 'active' : ''}" id="atom-btn">Atoms</button>
                </div>

                <!-- Spinner -->
                <div id="chsp-mode-spinner" class="chsp-mode-spinner spinner-border text-secondary" role="status" style="display: none;">
                    <span class="visually-hidden">Loading...</span>
                </div>
            </div>
        `;

        // Handle structure button
        const structureBtn = toggleContainer.querySelector('#structure-btn') as HTMLButtonElement;
        structureBtn.onclick = () => this._handleToggle('structure');

        // Handle atom button
        const atomBtn = toggleContainer.querySelector('#atom-btn') as HTMLButtonElement;
        atomBtn.onclick = () => this._handleToggle('atom');
        return toggleContainer.firstElementChild as HTMLElement;
    }

    /**
     * Handle toggle button click
     * @param mode flag indicating if the toggle should be checked
     */
    private _handleToggle(mode: DisplayMode): void {
        const isPerAtom = mode === 'atom';

        // Activate/desactivate structure button
        const structureBtn = this._toggleContainer.querySelector(
            '#structure-btn'
        ) as HTMLButtonElement;
        structureBtn.classList.toggle('active', !isPerAtom);

        // Activate/desactivate atom button
        const atomBtn = this._toggleContainer.querySelector('#atom-btn') as HTMLButtonElement;
        atomBtn.classList.toggle('active', isPerAtom);

        // Callback
        this.onchange(mode);
    }

    /**
     * Remove HTML added by DisplayToggle in the current document
     */
    remove(): void {
        this._shadow.host.remove();
    }

    /**
     * Toggle loading spinner
     * @param visible flag to toggle
     */
    loader(visible: boolean): void {
        // Show/hide spinnder
        const spinnerElement = this._toggleContainer.querySelector(
            '#chsp-mode-spinner'
        ) as HTMLDivElement;
        spinnerElement.style.display = visible ? 'inline-block' : 'none';

        // Toggle button visibility or disable state
        const buttons = this._toggleContainer.querySelectorAll('.btn-outline-secondary');
        buttons.forEach((button) => {
            if (visible) {
                button.setAttribute('disabled', 'true');
            } else {
                button.removeAttribute('disabled');
            }
        });
    }
}