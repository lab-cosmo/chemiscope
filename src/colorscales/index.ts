import {RGBColorMap} from './type';
import {Magma} from './magma';
import {Inferno} from './inferno';
import {Viridis} from './viridis';

import {ColorScale} from "../lib/plotly-scatter";

interface PlotlyColorMap {
    rgba: ColorScale;
    rgb: ColorScale;
}

function rgb_to_plotly(colormap: RGBColorMap): PlotlyColorMap {
    return {
        rgba: colormap.map((c) => {
            return [c[0], `rgba(${c[1][0]}, ${c[1][1]}, ${c[1][2]}, 0.5)`] as [number, string];
        }),
        rgb: colormap.map((c) => {
            return [c[0], `rgb(${c[1][0]}, ${c[1][1]}, ${c[1][2]})`] as [number, string];
        }),
    }
}

interface ColorMaps {
    [key: string]: PlotlyColorMap;
}

export const COLOR_MAPS: ColorMaps = {
    inferno: rgb_to_plotly(Inferno),
    magma: rgb_to_plotly(Magma),
    viridis: rgb_to_plotly(Viridis),
}
