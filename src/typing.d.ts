declare module '*.html' {
    const content: string;
    export default content;
}

declare module '*.svg' {
    const content: string;
    export default content;
}

declare module 'jsmol' {
    export interface JSmolApplet {
        _cover: (_: boolean) => void;
    }

    /**
     * Description of the Jmol Info object expected parameters as a typescript
     * interface.
     */
    export interface JmolInfo {
        use: 'HTML5' | 'Java' | 'WebGL';
        j2sPath: string;
        serverURL: string;
        height: number | string;
        width: number | string;
        color: string | [number, number, number];
        script: string;
        disableInitialConsole: boolean;
        disableJ2SLoadMonitor: boolean;
        zIndexBase: number;
        // other values are possible but not used in this project
    }

    /**
     * Declaration of the functions available on the Jmol global object
     */
    export interface JmolObject {
        setDocument(doc: boolean): void;
        getApplet(htmlId: string, info: Partial<JmolInfo>): JSmolApplet;
        script(applet: JSmolApplet, commands: string): void;
        getAppletHtml(applet: JSmolApplet): string;
        evaluateVar(applet: JSmolApplet, commands: string): unknown;
    }

    global {
        export interface Window {
            Jmol: JmolObject | undefined;
        }
    }
}
