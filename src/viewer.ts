import JSmolViewer from 'materials-cloud-viewer';

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

export class Viewer {
    private _structures: Structure[];
    private _viewer: JSmolViewer;

    constructor(id: string, j2sPath: string, structures: Structure[]) {
        this._viewer = new JSmolViewer(id, j2sPath);
        this._structures = structures;

        this.showStructure(0);
    }

    public showStructure(index: number) {
        if (index > this._structures.length) {
            console.warn(`invalid index in showStructure: got ${index}, max index is ${this._structures.length}`);
            return;
        }

        // TODO: does load remove the previous structures?
        this._viewer.load(`inline '${JSON.stringify(this._structures[index])}'`);
        const [a, b, c] = this._structures[index].cell;
        const a_str = `{${a[0]} ${a[1]} ${a[2]}}`;
        const b_str = `{${b[0]} ${b[1]} ${b[2]}}`;
        const c_str = `{${c[0]} ${c[1]} ${c[2]}}`;
        this._viewer.script(`unitcell [{0 0 0} ${a_str} ${b_str} ${c_str}];`);
    }
}
