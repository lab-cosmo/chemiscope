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
    cell: number[][],
}

export class JSmolViewer {
    private static _JSMOL_INFO = {
        use: "HTML5",
        width: 500,
        height: 500,
        script: "set antialiasdisplay; set frank off; set perspectiveDepth off",
        disableInitialConsole: true,
        disableJ2SLoadMonitor: true,
    };

    private _name: string;
    private _structures: Structure[];
    private _index: number;
    private _root: HTMLElement;
    private _slider: HTMLInputElement;
    private _slider_label: HTMLElement;
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

        const container = document.createElement("div");
        container.className = "sketchviz-viewer";
        this._root.appendChild(container);

        this._slider_label = document.createElement("span");
        this._slider_label.innerHTML = "Select an environment:"
        this._slider_label.style.width = "200px";
        this._slider_label.style.display = "inline-block";
        this._slider_label.style.margin = "20px";
        container.appendChild(this._slider_label);

        this._slider = document.createElement("input") as HTMLInputElement;
        this._slider.type = "range";
        this._slider.style.width = "300px";
        container.appendChild(this._slider);

        this._slider.oninput = () => {
            const i = this._slider.value;
            this._slider_label.innerHTML = `Environment ${i}`;

            // Default style
            this.script("select all;")
                this.script("wireframe 0.1; spacefill off; dots off;")
                this.script("color atoms translucent cpk;")

            // Style for the atoms in the environment (approximated as a 3.5A
            // cutoff)
            this.script(`select within(3.5, atomno = ${i});`)
                this.script("wireframe 0.13; spacefill 20%; dots off;")
                this.script("color atoms cpk;")

            // Style for the central atom
            this.script(`select atomno = ${i};`)
                this.script("wireframe off; spacefill 20%; dots 0.4;")
                this.script("color atoms green;")
        }

        const div = document.createElement("div");
        container.appendChild(div);

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
            return;
        }

        this._index = index;
        // TODO: does load remove the previous structures?
        this.script(`load inline '${JSON.stringify(this._structures[this._index])}';`);
        const [a, b, c] = this._structures[this._index].cell;
        const a_str = `{${a[0]} ${a[1]} ${a[2]}}`;
        const b_str = `{${b[0]} ${b[1]} ${b[2]}}`;
        const c_str = `{${c[0]} ${c[1]} ${c[2]}}`;
        this.script(`unitcell [{0 0 0} ${a_str} ${b_str} ${c_str}];`);

        this._slider_label.innerHTML = "Select an environment:"
        this._slider.min = "0";
        this._slider.max = `${this._structures[this._index].mol.a.length}`;
        this._slider.value = "0";
    }

    public script(commands: string) {
        this._Jmol.script(this._applet, commands);
    }
}
