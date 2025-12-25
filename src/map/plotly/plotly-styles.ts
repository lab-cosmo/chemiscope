import './plotly-scatter'; // Necessary for Plotly to add a <style> element
import type { PlotlyScatterElement } from './plotly-scatter';

declare global {
    interface CSSRule {
        selectorText: string;
    }
}

export function getPlotStyleSheet(plot: PlotlyScatterElement): CSSStyleSheet {
    return getStyleSheet(plot._fullLayout._modeBar._uid);
}

export const globalStyleSheet = getStyleSheet('global');

document.adoptedStyleSheets = [...document.adoptedStyleSheets, getDocumentStyleSheet()];

function getStyleSheet(name: string): CSSStyleSheet {
    const style = document.getElementById(`plotly.js-style-${name}`);

    // If the style element is missing (common in v3),
    // return an empty stylesheet instead of crashing.
    if (!style || !(style instanceof HTMLStyleElement) || !style.sheet) {
        return new CSSStyleSheet();
    }

    const sheet = new CSSStyleSheet();

    for (const rule of style.sheet.cssRules) {
        sheet.insertRule(rule.cssText);
    }

    style.remove();

    return sheet;
}

// Only keep the rules regarding .plotly-notifier which is a direct child of <body>.
function getDocumentStyleSheet(): CSSStyleSheet {
    const sheet = new CSSStyleSheet();

    for (const rule of globalStyleSheet.cssRules) {
        if (rule.selectorText.startsWith('.plotly-notifier')) {
            sheet.insertRule(rule.cssText);
        }
    }

    return sheet;
}
