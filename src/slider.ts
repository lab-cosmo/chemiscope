export class EnvironementSlider {
    private _root: HTMLElement;
    private _slider: HTMLInputElement;
    private _onchange: (index: number) => void;

    constructor(id: string, max: number) {
        const root = document.getElementById(id);
        if (root === null) {
            throw Error(`could not find HTML element #${id}`)
        }
        this._root = root;
        this._onchange = () => {}

        // TODO: Change to environnement when needed
        const centerType = 'structure';

        this._root.innerHTML = `
        <div class="input-group input-group-sm">
            <div class="input-group-prepend">
                <label class="input-group-text">${centerType}</label>
                <span class="input-group-text"><div class="skv-play-button"></div></span>
            </div>
            <input class="form-control custom-range" type='range' min=0 value=0 step=1></input>
        </div>
        `;

        this._slider = this._root.querySelector('input')!;
        this._slider.max = max.toString();
        this._slider.onchange = () => {
            this._onchange(parseInt(this._slider.value));
        }

        const play = this._root.querySelector('.skv-play-button')! as HTMLElement;
        const installTimeout = () => {
            setTimeout(() => {
                if (play.classList.contains('skv-pause')) {
                    // change value
                    let value = parseInt(this._slider.value);
                    value += 1;
                    if (value > max) {
                        value = 0;
                    }
                    this._slider.value = value.toString();
                    this._onchange(value);
                    // reinstall timeout
                    installTimeout();
                }
            }, 750)
        }

        play.onclick = () => {
            play.classList.toggle('skv-pause');
            installTimeout();
        }
    }

    /// Call the given callback when the user change the slider value
    public onChange(callback: (index: number) => void) {
        this._onchange = callback;
    }

    /// The environement index changed outside, update the slider
    public changed(index: number) {
        this._slider.value = index.toString();
    }
}
