/**
 * Draw each symbol using only the points that belong to it
 *
 * @packageDocumentation
 * @module map
 */

interface IndexBufferOptions {
    data: number[];
    type: 'uint32';
    primitive: 'points';
}

interface IndexBuffer {
    destroy: () => void;
}

/** One entry per point, nonzero when the point uses this symbol */
interface SymbolMask {
    data: Uint8Array;
    destroy: () => void;
}

interface ScatterRenderer {
    regl: { elements: (options: IndexBufferOptions) => IndexBuffer };
    markerTextures: unknown[];
}

interface PointGroup {
    // spatial index used by the renderer drawing path
    tree?: unknown;
    // one mask per symbol or true if every point uses that symbol
    activation: (SymbolMask | boolean | null | undefined)[];
}

interface SymbolPoints {
    indexBuffer: IndexBuffer;
    pointCount: number;
}

// weak keys let old masks be collected after the renderer replaces them
const SYMBOL_BUFFERS = new WeakMap<SymbolMask, SymbolPoints>();

/** Build drawing options that submit only this symbol's points */
export function getSymbolDrawOptions(
    renderer: ScatterRenderer,
    pointGroup: PointGroup,
    symbolId: number,
    selectedPointIds: unknown
): object | null {
    const symbolMask = pointGroup.activation[symbolId];

    // keep the original handling for clusters, selections and uniform symbols
    if (pointGroup.tree || selectedPointIds || !symbolMask || symbolMask === true) {
        return null;
    }

    // reuse the index buffer while the symbol mask is unchanged
    let symbolPoints = SYMBOL_BUFFERS.get(symbolMask);
    if (symbolPoints === undefined) {
        // turn a mask such as [0, 1, 0, 1] into point ids [1, 3]
        const pointIds: number[] = [];
        for (let i = 0; i < symbolMask.data.length; i++) {
            if (symbolMask.data[i] !== 0) {
                pointIds.push(i);
            }
        }

        // point ids can exceed 65535 when lod is off
        const indexBuffer = renderer.regl.elements({
            data: pointIds,
            type: 'uint32',
            primitive: 'points',
        });
        symbolPoints = { indexBuffer, pointCount: pointIds.length };
        SYMBOL_BUFFERS.set(symbolMask, symbolPoints);

        // the renderer destroys the mask when symbols change or the plot is removed
        // release our gpu buffer too, garbage collection cannot free it
        const destroyMask = symbolMask.destroy;
        symbolMask.destroy = () => {
            indexBuffer.destroy();
            destroyMask();
        };
    }

    return {
        // keep the existing positions, colors, sizes and other drawing settings
        ...pointGroup,
        markerTexture: renderer.markerTextures[symbolId],

        // every submitted point now belongs to this symbol, so no mask is needed
        activation: true,
        elements: symbolPoints.indexBuffer,
        count: symbolPoints.pointCount,
        // start at the beginning of this symbol's point-id list
        offset: 0,
    };
}
