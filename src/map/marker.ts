import assert from 'assert';

import { Indexes } from '../indexer';

/// Data associated with markers on the map indicating selected environments/structures
export class MarkerData {
    /// color of the marker, synchronized with the ViewersGrid where appropriate
    public color: string;
    /// index of the currently displayed structure
    public current: number;
    /// Marker indicating the position of the point in 2D mode
    /// Using an HTML element is much faster than using Plolty to restyle the full plot,
    /// especially with more than 100k points. In 3D mode, a separate trace is
    /// used instead.
    public marker: HTMLElement;

    constructor(guid: string, color: string, environment: number, visible: boolean = true) {
        const marker = document.createElement('div');
        marker.classList.add('chsp-structure-marker');
        marker.id = `chsp-selected-${guid}`;
        marker.style.backgroundColor = color;
        this.marker = marker;

        this.toggleVisible(visible);

        this.color = color;
        this.current = environment;
    }

    public update(x: number, y: number): void {
        if (!isFinite(x) || !isFinite(y)) {
            this.toggleVisible(false);
            return;
        }
        this.marker.style.top = `${y}px`;
        this.marker.style.right = `${x}px`;
    }
    public select(indexes: Indexes): boolean {
        if (this.current !== indexes.environment) {
            this.current = indexes.environment;
            return true;
        } else {return false; }
    }
    public activate(): void {
        this.marker.classList.toggle('chsp-active-structure', true);
    }
    public deactivate(): void {
        this.marker.classList.toggle('chsp-active-structure', false);
    }
    public remove(): void {
        assert(this.marker.parentNode !== null);
        this.marker.parentNode.removeChild(this.marker);
    }

    public toggleVisible(visible: boolean = false): void {
        if (visible) {
            this.marker.style.display = 'block';
        } else {
            this.marker.style.display = 'none';
        }
    }
}
