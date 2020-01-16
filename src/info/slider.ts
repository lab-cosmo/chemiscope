import {Target} from '../dataset';

/// A simple slider to select an environment / play the trajectory
export class Slider {
    public onchange: () => void;
    private _slider: HTMLInputElement;
    private _play: HTMLElement;
    private _delay: HTMLInputElement;

    constructor(root: HTMLElement, target: Target, delay: HTMLInputElement) {
        const template = document.createElement('template');
        template.innerHTML = `<div class="input-group input-group-sm">
            <div class="input-group-prepend">
                <span class="input-group-text"><div class="skv-play-button"></div></span>
            </div>
            <input class="form-control custom-range skv-${target}-range" type='range' min=0 value=0 step=1></input>
        </div>`;
        const group = template.content.firstChild! as HTMLElement;
        root.appendChild(group);

        this._slider = group.querySelector('input')! as HTMLInputElement;
        this._play = group.querySelector('.skv-play-button')! as HTMLElement;
        this._delay = delay;

        this._play.onclick = () => {
            this._play.classList.toggle('skv-playing');
            this._step();
        }

        this._slider.onchange = () => this.onchange();
        this.onchange = () => {};
    }

    public reset(max: number) {
        this._slider.value = "0";
        this._slider.max = max.toString();
    }

    public update(value: number) {
        this._slider.value = value.toString();
    }

    public value(): number {
        return parseInt(this._slider.value);
    }

    private _step() {
        setTimeout(() => {
            if (this._play.classList.contains('skv-playing')) {
                const value = (this.value() + 1) % parseInt(this._slider.max);
                this.update(value);
                this.onchange();
                // contibue playing until the 'skv-playing' class
                // is no longer there
                this._step();
            }
        }, parseInt(this._delay.value) * 100);
    }
}
