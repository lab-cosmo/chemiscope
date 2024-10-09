/**
 * Component for managing modals, with two main roles: apply classes (and
 * eventually add/remove from the dom) to the modal and its backdrop at the
 * right time, and listen for events causing the modal to close (e.g. clicking
 * on the backdrop).
 *
 * The classes and styles applied match those used by Bootstrap. They are
 * summarized in this table:
 *               | Modal         | Backdrop                  |
 *  Closed       | display: none | <removed>                 |
 *  Closing      | -             | .modal-backdrop.fade      |
 *  Open/Opening | .show         | .modal-backdrop.fade.show |
 * More classes must be set by the user for correct styling, such as .modal for
 * the modal. See [the Bootstrap docs](https://getbootstrap.com/docs/5.2/components/modal/)
 * for details.
 *
 * The modal is inserted in the body and is contained inside a shadow
 * root. It can be styled as follows:
 *  modal.shadow.adoptedStyleSheets = [...]
 *
 * When the modal is closed, the focus, if any, is restored to the element that
 * was focused prior to opening the modal.
 */
export default class Modal {
    private _activeElement: Element | null = null;
    private _backdrop: HTMLElement;
    private _open = false;

    public shadow: ShadowRoot;

    /**
     * Create a Modal instance.
     * @param _element The root element to use as a modal. It should have the
     *  .modal class.
     */
    constructor(private _element: HTMLElement) {
        // Create a shadow root in the body to hold the modal.
        const hostElement = document.createElement('div');
        document.body.appendChild(hostElement);

        this.shadow = hostElement.attachShadow({ mode: 'open' });
        this.shadow.appendChild(this._element);

        // Initialize the backdrop element. It will be inserted and removed from
        // the DOM following the modal's lifecycle.
        this._backdrop = document.createElement('div');
        this._backdrop.classList.add('modal-backdrop', 'fade');

        // Listen for click events on the root element (excluding its children),
        // which corresponds to the backdrop.
        this._element.addEventListener('click', (event) => {
            if (event.target === event.currentTarget) {
                this.close();
            }
        });

        // Find close buttons inside the modal.
        for (const el of this._element.querySelectorAll('[data-bs-dismiss="modal"]')) {
            el.addEventListener('click', () => {
                this.close();
            });
        }

        // Listen for Escape key presses.
        document.addEventListener('keydown', (event) => {
            if (event.key === 'Escape') {
                event.preventDefault();
                this.close();
            }
        });

        // Listeners called once the modal is closed

        this._element.addEventListener('transitionend', (event) => {
            if (event.target === event.currentTarget && !this._open) {
                this._element.style.setProperty('display', 'none');
            }
        });

        this._backdrop.addEventListener('transitionend', () => {
            if (!this._open) {
                this._backdrop.remove();

                if (this._activeElement && this._activeElement instanceof HTMLElement) {
                    this._activeElement.focus();
                }

                this._activeElement = null;
            }
        });
    }

    close(): void {
        this._backdrop.classList.remove('show');
        this._element.classList.remove('show');
        this._open = false;
    }

    open(): void {
        this._element.getRootNode().appendChild(this._backdrop);
        this._element.style.setProperty('display', 'block');

        // Trick to allow transitioning immediately after setting display to block
        // eslint-disable-next-line @typescript-eslint/no-unused-vars
        const _ = this._element.offsetHeight;

        this._backdrop.classList.add('show');
        this._element.classList.add('show');

        // Unfocus the active element and store it in order to focus it again
        // once the modal is closed.
        this._activeElement = getActiveElement();

        if (this._activeElement && this._activeElement instanceof HTMLElement) {
            this._activeElement.blur();
        }

        this._open = true;
    }

    toggle(value: boolean = !this.open): void {
        if (value) {
            this.open();
        } else {
            this.close();
        }
    }

    remove(): void {
        this.shadow.host.remove();
    }
}

// Returns the active element in `root`, if any, no matter how nested it is in
// shadow roots. When reaching a closed shadow root containing that element,
// the function will return the root's host element.
function getActiveElement(root: Document | ShadowRoot = document): Element | null {
    const activeElement = root.activeElement;

    return activeElement?.shadowRoot ? getActiveElement(activeElement.shadowRoot) : activeElement;
}
