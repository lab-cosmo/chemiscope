import assert from 'assert';

import {EnvironmentIndexer, Indexes} from '../indexer';

/// interface to contain synchronized parameters for marker instances
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

    constructor(
        guid: string,
        color: string,
        environment: number,
        visible: boolean = true,
    ) {
      const marker = document.createElement('div');
      marker.classList.add('chsp-structure-marker');
      marker.id = `chsp-selected-${guid}`;
      marker.style.backgroundColor = color;
      this.marker = marker;

      this.toggleVisible(visible);

      this.color = color;
      this.current = environment;
    }

    public update(x: number, y: number, bounds: number[][]) {
          if (!isFinite(x) || !isFinite(y)) {
              this.toggleVisible(false);
          } else if (x < bounds[0][0] || x > bounds[0][1] || y < bounds[1][0] || y > bounds[1][1]) {
              this.toggleVisible(false);
          } else {
              this.toggleVisible(true);
          }
          this.marker.style.top = `${y}px`;
          // const plotWidth = this._plot.getBoundingClientRect().width;
          const plotWidth = bounds[0][1] - bounds[0][0];
          this.marker.style.right = `${plotWidth - x}px`;

    }
    public select(indexes: Indexes): void {
      if (this.current !== indexes.environment) {
          this.current = indexes.environment;
          // this.update();
      }
    }
    public activate() {
      this.marker.classList.toggle('chsp-active-structure', true);
    }
    public deactivate() {
      this.marker.classList.toggle('chsp-active-structure', false);
    }
    public remove() {
      assert(this.marker.parentNode !== null);
      this.marker.parentNode.removeChild(this.marker);
    }

    public toggleVisible(visible: boolean = false) {
      if (visible) {
          this.marker.style.display = 'block';
      } else {
          this.marker.style.display = 'none';
      }
    }

}
