import {EnvironmentIndexer, Indexes} from '../indexer';

/// interface to contain synchronized parameters for marker instances
class MarkerData {
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
        current: number,
        marker: HTMLElement,
        color?: string,
    ) {
      this.color = color;
      this.current = current;
      this.marker = marker;
    }

    public update() {}
    public select(indices: Indexes) {}
    public activate() {}
    public deactivate() {}
    public remove() {}

}
