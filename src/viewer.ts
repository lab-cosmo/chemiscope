interface Atom {
    l: string;
    x: number;
    y: number;
    z: number;
}

interface Structure {
    mol: {
        a: Atom[],
    };
}

export class JSmolViewer {
    private static _JSMOL_INFO = {
        use: "HTML5",
        width: 500,
        height: 500,
        script: "set antialiasdisplay; set frank off",
        disableInitialConsole: true,
        disableJ2SLoadMonitor: true,
    };

    private _name: string;
    private _structures: Structure[];
    private _root: HTMLElement;
    private _applet: any;
    private _Jmol: any;

    constructor(name: string, structures: Structure[]) {
        if (!("Jmol" in window)) {
            console.error("Jmol is required. Please load it.");
        }
        this._Jmol = window["Jmol"];
        this._name = name;
        this._structures = structures;
    }

    public setup(root: HTMLElement, j2sPath: string) {
        this._root = root;

        const div = document.createElement("div");
        this._root.appendChild(div);

        this._applet = this._Jmol.getApplet(this._name + "-JmolApplet", {
            ...JSmolViewer._JSMOL_INFO,
            j2sPath: j2sPath,
        });

        div.innerHTML = this._Jmol.getAppletHtml(this._applet);
        // Jmol rely on this script being implicitly executed, but this is not
        // the case when using innerHTML (compared to jquery .html()). So let's
        // manually execute it
        this._applet._cover(false);

        this.showStructure(0);
    }

    public showStructure(index: number) {
        if (index > this._structures.length) {
            console.warn(`invalid index in showStructure: got ${index}, max index is ${this._structures.length}`);
        }
        // TODO: does load remove the previous structures?
        this.script(`load inline '${JSON.stringify(this._structures[index])}'`);
        this.script("spin on");
    }

    public script(commands: string) {
        this._Jmol.script(this._applet, commands);
    }
}
