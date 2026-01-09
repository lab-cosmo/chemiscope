/**
 * @packageDocumentation
 * @module info
 */

import assert from 'assert';

import { Target } from '../dataset';
import { binarySearch } from '../utils';

/** @hidden
 * A simple slider to select an environment / play the trajectory
 */
export class Slider {
    /** Callback fired when the use changes the slider value */
    public onchange: () => void;
    /**
     * Callback fired when the use click the play button. The `advance`
     * callback indicate whether to continue playback or not.
     */
    public startPlayback: (advance: () => boolean) => void;

    private _slider: HTMLInputElement;
    private _play: HTMLElement;
    private _valid: number[];

    /**
     * Create and append a new slider inside the given `HTMLElement`.
     *
     * @param root   where to append the new slider
     * @param target is this slider related to atom or structure
     */
    constructor(root: HTMLElement, target: Target) {
        const template = document.createElement('template');
        template.innerHTML = `<div part="chsp-slider" class="input-group input-group-sm">
            <span class="input-group-text"><div class="chsp-play-button"></div></span>
            <input class="form-control form-range chsp-${target}-slider" type='range' min=0 value=0 step=1 style="height: 1.9rem">
        </div>`;
        const group = template.content.firstChild as HTMLElement;
        root.appendChild(group);

        this._valid = [];
        this._slider = group.querySelector('input') as HTMLInputElement;
        this._play = group.querySelector('.chsp-play-button') as HTMLElement;

        this._play.onclick = () => {
            this._play.classList.toggle('chsp-playing');
            this.startPlayback(() => {
                return this._play.classList.contains('chsp-playing');
            });
        };

        this._slider.onchange = () => this.onchange();
        this.onchange = () => {};
        this.startPlayback = () => {};
    }

    /**
     * Reset the slider to the first entry in the given `valid` values
     */
    public reset(valid: number[]): void {
        this._valid = valid;
        this._slider.value = '0';
        this._slider.max = `${this._valid.length - 1}`;
    }

    /**
     * Update the slider, changing its value to the given one.
     */
    public update(value: number): void {
        const position = binarySearch(this._valid, value);
        assert(position !== -1);
        this._slider.value = position.toString();
    }

    /**
     * Get the current value of the slider
     */
    public value(): number {
        return this._valid[parseInt(this._slider.value, 10)];
    }

    /**
     * Is the slider currently playing?
     */
    public get isPlaying(): boolean {
        return this._play.classList.contains('chsp-playing');
    }
}
