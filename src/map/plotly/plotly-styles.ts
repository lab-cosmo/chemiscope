import assert from 'assert';

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
  assert(style instanceof HTMLStyleElement);
  assert(style.sheet);


  const sheet = new CSSStyleSheet();

  for (const rule of Array.from(style.sheet.cssRules)) {
    sheet.insertRule(rule.cssText);
  }

  style.remove();

  return sheet;
}

// Only keep the rules regarding .plotly-notifier which is a direct child of <body>.
function getDocumentStyleSheet(): CSSStyleSheet {
    const sheet = new CSSStyleSheet();

    for (const rule of Array.from(globalStyleSheet.cssRules)) {
        if (rule.selectorText.startsWith('.plotly-notifier')) {
            sheet.insertRule(rule.cssText);
        }
    }

    return sheet;
}
