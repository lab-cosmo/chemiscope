import { Streamlit } from 'streamlit-component-lib';

const ROOT_ID = 'chemiscope-root';

export const CONFIG = {
    map: 'chemiscope-map',
    info: 'chemiscope-info',
    meta: 'chemiscope-meta',
    structure: 'chemiscope-structure',
} as const;

export function getOrCreateRoot(): HTMLDivElement {
    let root = document.getElementById(ROOT_ID) as HTMLDivElement | null;
    if (!root) {
        root = document.createElement('div');
        root.id = ROOT_ID;
        root.style.width = '100%';
        root.style.display = 'flex';
        document.body.appendChild(root);
    }
    return root;
}

export function applyWidthPolicy(widthArg: string | number, root: HTMLDivElement): void {
    const viewer = root.querySelector('.chemiscope-streamlit') as HTMLElement;
    if (!viewer) {
        return;
    }

    if (widthArg === 'stretch') {
        root.style.justifyContent = 'stretch';
        root.style.width = '100%';
        viewer.style.width = '100%';
    } else if (typeof widthArg === 'number') {
        root.style.width = '100%';
        root.style.justifyContent = 'center';
        viewer.style.width = widthArg + 'px';
    }
}

export function applyHeightPolicy(heightArg: number, root: HTMLElement): void {
    root.style.height = heightArg + 'px';
    Streamlit.setFrameHeight(heightArg);
}

export function toggleLoadingVisible(visible: boolean = true): void {
    const loader = document.getElementById(`${ROOT_ID}-loading`);
    if (loader) {
        loader.style.display = visible ? 'flex' : 'none';
    }
}

export function displayWarning(message: any, timeout: number = 4000): void {
    const display = document.getElementById(`${ROOT_ID}-warning-display`);
    if (!display) {
        console.error('Warning display element not found:', `${ROOT_ID}-warning-display`);
        return;
    }

    const textMessage = message instanceof Error ? message.toString() : String(message);
    const paragraph = display.querySelector('p');
    if (paragraph) {
        paragraph.innerText = textMessage;
    }

    display.style.display = 'flex';

    if (timeout > 0) {
        setTimeout(() => {
            display.style.display = 'none';
        }, timeout);
    }
}

export function generateHTMLForMode(mode: string, noInfo: boolean = false): string {
    let layout: string;

    const mapStyle = mode === 'structure' ? 'display: none;' : '';
    const structureStyle = mode === 'map' ? 'display: none;' : '';
    const infoStlye = noInfo ? 'display: none;' : '';

    if (mode === 'structure') {
        layout = `
            <div class="visualizer-container visualizer-structure-mode">
                <div class="visualizer-column">
                    <div id="${CONFIG.meta}"></div>
                    <div id="${CONFIG.structure}" class="visualizer-item"></div>
                    <div id="${CONFIG.info}" class="visualizer-info" style="${infoStlye}"></div>
                    <div id="${CONFIG.map}" style="${mapStyle}"></div> </div>
            </div>`;
    } else if (mode === 'map') {
        layout = `
            <div class="visualizer-container visualizer-map-mode">
                <div class="visualizer-column">
                    <div id="${CONFIG.meta}"></div>
                    <div id="${CONFIG.map}" class="visualizer-item"></div>
                    <div id="${CONFIG.info}" class="visualizer-info" style="${infoStlye}"></div>
                    <div id="${CONFIG.structure}" style="${structureStyle}"></div> </div>
            </div>`;
    } else {
        layout = `
            <div class="visualizer-container">
                <div class="visualizer-column-right">
                    <div id="${CONFIG.meta}"></div>
                    <div id="${CONFIG.map}" class="visualizer-item"></div>
                </div>
                <div class="visualizer-column">
                    <div id="${CONFIG.structure}" class="visualizer-item"></div>
                    <div id="${CONFIG.info}" class="visualizer-info" style="${infoStlye}"></div>
                </div>
            </div>`;
    }

    return `<div class="chemiscope-streamlit">${layout}</div>`;
}
