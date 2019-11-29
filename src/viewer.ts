import JSmolViewer from 'materials-cloud-viewer';

export class StructureViewer {
    private _structures: string[];
    private _viewer: JSmolViewer;

    constructor(id: string, j2sPath: string, structures: string[]) {
        this._viewer = new JSmolViewer(id, j2sPath);
        this._structures = structures;

        this.showStructure(0);
    }

    public changeDataset(structures: string[]) {
        this._structures = structures;
        this.showStructure(0);
    }

    public showStructure(index: number) {
        if (index > this._structures.length) {
            console.warn(`invalid index in showStructure: got ${index}, max index is ${this._structures.length}`);
            return;
        }

        this._viewer.load(`inline '${this._structures[index]}'`);
    }
}
