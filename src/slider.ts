/// A simple slider to select an environment / play the trajectory
export class Slider {
    public root: HTMLElement;
    public onchange: () => void;
    private _name: string;
    private _slider: HTMLInputElement;
    private _play: HTMLElement;
    private _label: HTMLLabelElement;

    constructor(name: string) {
        const template = document.createElement('template');
        template.innerHTML = `<div class="input-group input-group-sm">
            <div class="input-group-prepend">
                <label class="input-group-text skv-slider-label">${name}</label>
                <span class="input-group-text"><div class="skv-play-button"></div></span>
            </div>
            <input class="form-control custom-range" type='range' min=0 value=0 step=1></input>
        </div>`;
        this._name = name;
        this.root = template.content.firstChild! as HTMLElement;
        this._slider = this.root.querySelector('input')! as HTMLInputElement;
        this._play = this.root.querySelector('.skv-play-button')! as HTMLElement;
        this._label = this.root.querySelector('label')! as HTMLLabelElement;

        this._play.onclick = () => {
            this._play.classList.toggle('skv-playing');
            this.step();
        }

        this._slider.onchange = () => this.onchange();
        this.onchange = () => {};
    }

    public reset(max: number) {
        this._label.innerText = `${this._name} 0`;
        this._slider.value = "0";
        this._slider.max = max.toString();
    }

    public update(value: number) {
        this._label.innerText = `${this._name} ${value}`;
        this._slider.value = value.toString();
    }

    public value(): number {
        return parseInt(this._slider.value);
    }

    public step() {
        setTimeout(() => {
            if (this._play.classList.contains('skv-playing')) {
                const value = (this.value() + 1) % parseInt(this._slider.max);
                this.update(value);
                this.onchange();
                // contibue playing until the 'skv-playing' class
                // is no longer there
                this.step();
            }
        }, 750)
    }
}
