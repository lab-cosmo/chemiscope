/**
 * @packageDocumentation
 * @module info
 */

import {Target} from '../dataset';

/** @hidden
 * A simple slider to select an environment / play the trajectory
 */
export class Slider {
    private _slider: HTMLInputElement;
    private _play: HTMLElement;

    /** Callback fired when the use changes the slider value */
    public onchange: () => void;
    /**
     * Callback fired when the use click the play button. The `advance`
     * callback indicate whether to continue playback or not.
     */
    public startPlayback: (advance: () => boolean) => void;

    /**
     * Create and append a new slider inside the given `HTMLElement`.
     *
     * @param root   where to append the new slider
     * @param target is this slider related to atom or structure
     */
    constructor(root: HTMLElement, target: Target) {
        const template = document.createElement('template');
        template.innerHTML = `<div class="input-group input-group-sm">
            <div class="input-group-prepend">
                <span class="input-group-text"><div class="chsp-play-button"></div></span>
            </div>
            <input class="form-control custom-range chsp-${target}-slider" type='range' min=0 value=0 step=1></input>
        </div>`;
        const group = template.content.firstChild! as HTMLElement;
        root.appendChild(group);

        this._slider = group.querySelector('input')! as HTMLInputElement;
        this._play = group.querySelector('.chsp-play-button')! as HTMLElement;

        this._play.onclick = () => {
            this._play.classList.toggle('chsp-playing');
            this.startPlayback(() => {return this._play.classList.contains('chsp-playing');});
        }

        this._slider.onchange = () => this.onchange();
        this.onchange = () => {};
        this.startPlayback = () => {};
    }

    /**
     * Reset the slider to 0, and update the maximal value to `max`
     */
    public reset(max: number) {
        this._slider.value = "0";
        this._slider.max = max.toString();
    }

    /**
     * Update the slider, changing its value to the given one.
     */
    public update(value: number) {
        this._slider.value = value.toString();
    }

    /**
     * Get the current value of the slider
     */
    public value(): number {
        return parseInt(this._slider.value);
    }
}
