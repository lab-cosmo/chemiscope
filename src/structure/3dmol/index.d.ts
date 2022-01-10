interface ElementColors {
    [name: string]: ColorSpec;
}

export interface ViewerSpec {
    /** Callback function to be immediately executed on this viewer */
    callback: () => void; // TODO: actual type
    /** Object defining default atom colors as atom => color property value
     * pairs for all models within this viewer */
    defaultcolors: ElementColors;
    /** Whether to disable mouse handlers */
    nomouse: boolean;
    /** Color of the canvas' background */
    backgroundColor: ColorSpec;
    backgroundAlpha: number;

    disableFog: boolean;

    camerax: number;
    hoverDuration: number;

    /** id of the canvas */
    id: string;
    /** default 5 */
    cartoonQuality: number;

    row: number;
    col: number;
    rows: number;
    cols: number;
    canvas: unknown;
    viewers: unknown;

    minimumZoomToDistance: unknown;
    lowerZoomLimit: unknown;
    upperZoomLimit: unknown;

    antialias: boolean;
    control_all: boolean;
    orthographic: boolean;
}

export interface ParserOptionsSpec {
    /** true if you want to add to a new frame and false otherwise; supported by all */
    frames: boolean;
    /** object specifying the vibration behavior; supported by all */
    vibrate: {
        /** number of frames to be created, default to 10; supported by all */
        frames: number;
        /** amplitude of distortion, default to 1 (full); supported by all */
        amplitude: number;
    };
    /** specifies weather or not multiple models are being defined; supported by xyz, sdf, or mol2 */
    multimodel: boolean;
    /** specifies weather or not the model is of one molecule; supported by xyz, sdf, mol2 */
    onemol: boolean;
    /** do not strip hydrogens; supported by sdf, mol2 */
    keepH: boolean;
    /** used to define ChemDoodle styles; supported by cdjson */
    parseStyle: unknown;
    /** dictating weather or not to do assembly; supported by mcif */
    doAssembly: boolean;
    /** Set to true if you wish to duplicate assembly atoms otherwise false;
     * supported by all formats with symmetries. Not duplicating will result
     * in faster rendering but it will not be possible to individually style
     * symmetries. */
    duplicateAssemblyAtoms: boolean;
    /** shift symmetry mates so their centroid is in the unit cell */
    normalizeAssembly: boolean;
    /** do not detect bonds between symmetries generated with
     * duplicateAssemblyAtoms (cif only - other formats never make bonds
     * between symmetries) */
    dontConnectDuplicatedAtoms: boolean;
    /** boolean dictating the presence of a secondary structure; supported by pdb */
    noSecondaryStructure: boolean;
    /** do not compute ss; supported by pdb */
    noComputeSecondaryStructure: boolean;
    /** which alternate location to select, if present; '*' to load all;
     * supported by pdb */
    altLoc: string;
    /** index of the assembly in symmetry; supported by mmtf */
    assemblyIndex: number;
}

export type ColorSpec = string | number;

export interface AtomSpec {
    index: number;
    /** Parent residue name */
    resn: string;
    /** Atom's x coordinate */
    x: number;
    /** Atom's y coordinate */
    y: number;
    /** Atom's z coordinate */
    z: number;
    /** Atom's color, as hex code or built-in color string */
    color: ColorSpec;
    /** Hex code for color to be used for surface patch over this atom */
    surfaceColor: ColorSpec;
    /** Element abbreviation (e.g. 'H', 'Ca', etc) */
    elem: string;
    /** Set to true if atom is a heteroatom */
    hetflag: boolean;
    /** Chain this atom belongs to, if specified in input file (e.g 'A' for chain A) */
    chain: string;
    /** Residue number */
    resi: number;
    icode: number;
    rescode: number;
    /** Atom's serial id number */
    serial: number;
    /** Atom name; may be more specific than 'elem' (e.g 'CA' for alpha carbon) */
    atom: string;
    /** Array of atom ids this atom is bonded to */
    bonds: number[];
    /** Secondary structure identifier (for cartoon render; e.g. 'h' for helix) */
    ss: string;
    /** true if this atom forms only single bonds or no bonds at all */
    singleBonds: boolean;
    /** Array of this atom's bond orders, corresponding to bonds identfied by 'bonds' */
    bondOrder: number[];
    /** Optional mapping of additional properties */
    properties: unknown;
    /** Atom b factor data */
    b: number;
    /** If applicable, this atom's record entry from the input PDB file (used to output new PDB from models) */
    pdbline: string;
    /** Set this flag to true to enable click selection handling for this atom */
    clickable: boolean;
    // TODO
    // /** Callback click handler function to be executed on this atom and its parent viewer */
    // callback: function;
    /** for selection, inverts the meaning of the selection */
    invert: boolean;
    altLoc: string;
}

export interface AtomSelectionSpec extends Omit<AtomSpec, 'bonds'> {
    /** a single model or list of models from which atoms should be selected. Can also specify by numerical creation order. Reverse indexing is allowed (-1 specifies last added model). */
    model: GLModel | number;
    /** overloaded to select number of bonds, e.g. {bonds: 0} will select all nonbonded atoms */
    bonds: number;
    /** user supplied function that gets passed an {AtomSpec} and should return true if the atom should be selected */
    predicate: (spec: AtomSpec) => boolean;
    /** if set, inverts the meaning of the selection */
    invert: boolean;
    /** if set, expands the selection to include all atoms of any residue that has any atom selected */
    byres: boolean;
    /** expands the selection to include all atoms within a given distance from the selection */
    expand: number;
    /** intersects the selection with the set of atoms within a given distance from another selection */
    within: Partial<WithinSelectionSpec>;
    /** take the intersection of the provided lists of {AtomSelectionSpec}s */
    and: Array<Partial<AtomSelectionSpec>>;
    /** take the union of the provided lists of {AtomSelectionSpec}s */
    or: Array<Partial<AtomSelectionSpec>>;
    /** take the inverse of the provided {AtomSelectionSpec} */
    not: Partial<AtomSelectionSpec>;
}

export interface WithinSelectionSpec {
    /** the distance in angstroms away from the atom selection to include atoms in the parent selection */
    distance: number;
    /** if set, selects atoms not within distance range for intersection */
    invert: boolean;
    /** the selection of atoms against which to measure the distance from the parent atom selection */
    sel: Partial<AtomSelectionSpec>;
}

export interface AtomStyleSpec {
    /** draw bonds as lines */
    line: Partial<LineStyleSpec>;
    /** draw atoms as crossed lines (aka stars) */
    cross: Partial<CrossStyleSpec>;
    /** draw bonds as capped cylinders */
    stick: Partial<StickStyleSpec>;
    /** draw atoms as spheres */
    sphere: Partial<SphereStyleSpec>;
    /** draw cartoon representation of secondary structure */
    cartoon: Partial<CartoonStyleSpec>;
    /** invisible style for click handling only */
    clicksphere: Partial<ClickSphereStyleSpec>;
}

export interface LineStyleSpec {
    /** do not show line */
    hidden: boolean;
    /** deprecated due to vanishing browser support */
    linewidth: number;
    /** element based coloring */
    colorscheme: Partial<ColorschemeSpec>;
    /** fixed coloring, overrides colorscheme */
    color: ColorSpec;
    opacity: number;
}

export interface CrossStyleSpec {
    /** do not show */
    hidden: boolean;
    /** deprecated due to vanishing browser support */
    linewidth: number;
    radius: number;
    /** scale radius by specified amount */
    scale: number;
    /** element based coloring */
    colorscheme: Partial<ColorschemeSpec>;
    /** fixed coloring, overrides colorscheme */
    color: ColorSpec;
    opacity: number;
}

export interface StickStyleSpec {
    /** do not show */
    hidden: boolean;
    radius: number;
    /** draw all bonds as single bonds if set */
    singleBonds: boolean;
    /** element based coloring */
    colorscheme: Partial<ColorschemeSpec>;
    /** fixed coloring, overrides colorscheme */
    color: ColorSpec;
    opacity: number;
}

export interface SphereStyleSpec {
    /** do not show atom */
    hidden: boolean;
    /** override van der waals radius */
    radius: number;
    /** scale radius by specified amount */
    scale: number;
    /** element based coloring */
    colorscheme: Partial<ColorschemeSpec>;
    /** fixed coloring, overrides colorscheme */
    color: ColorSpec;
    opacity: number;
}

export interface ColorschemeSpec {
    // TODO missing property name in docs
    // color>Carbon - use default element colors but with carbon set to specify html color string
    /** PyMol secondary colorscheme */
    ssPyMOL: string;
    /** Jmol secondary colorscheme */
    ssJmol: string;
    /** Jmol primary colorscheme */
    Jmol: string;
    /** default colorscheme */
    default: string;
    /** amino acid colorscheme */
    amino: string;
    /** shapely protien colorscheme */
    shapely: string;
    /** nucleic acid colorscheme */
    nucleic: string;
    /** standard chain colorscheme */
    chain: string;
    /** chain Hetatm colorscheme */
    chainHetatm: string;
    /** atomSpec property. Example 'b'. See AtomSpec. */
    prop: string;
    // TODO
    // /** Allows the user to provide a gradient to the colorscheme. */
    // gradient: Gradient;
    // /** map of a certain AtomSpec property to a color. */
    // map: object;
    /** Allows the user to provide a mapping of elements to colors to the
     * colorscheme. This can be done with any properties, and not just 'elem'. */
    elementMap: ElementColors;
    // /** Allows the user to provide a function for setting the colorschemes. */
    // colorfunc: function;
}

export interface CartoonStyleSpec {
    /** strand color, may specify as 'spectrum' which will apply reversed gradient based on residue number */
    color: ColorSpec;
    /** style of cartoon rendering (trace, oval, rectangle (default), parabola, edged) */
    style: string;
    /** whether to use constant strand width, disregarding secondary structure; use thickness to adjust radius */
    ribbon: boolean;
    /** whether to add arrows showing beta-sheet directionality; does not apply to trace or ribbon */
    arrows: boolean;
    /** whether to display alpha helices as simple cylinders; does not apply to trace */
    tubes: boolean;
    /** cartoon strand thickness, default is 0.4 */
    thickness: number;
    /** cartoon strand width, default is secondary structure-dependent; does not apply to trace or ribbon */
    width: number;
    /** set opacity from 0-1; transparency is set per-chain
     * with a warning outputted in the event of ambiguity
     *
     * In nucleic acids, the base cylinders obtain their color from the
     * atom to which the cylinder is drawn, which is 'N1' for purines (resn:
     * 'A', 'G', 'DA', 'DG') and 'N3' for pyrimidines (resn: 'C', 'U', 'DC',
     * 'DT'). The different nucleobases can therefore be distinguished as
     * follows:*/
    opacity: number;
}

export interface UnitCellStyleSpec {
    /** line style used to draw box */
    box: Partial<LineStyleSpec>;
    /** arrow specification of "a" axis */
    astyle: Partial<ArrowSpec>;
    /** arrow specification of "b" axis */
    bstyle: Partial<ArrowSpec>;
    /** arrow specification of "c" axis */
    cstyle: Partial<ArrowSpec>;
    /** label for a axis */
    alabel: string;
    /** label style for a axis */
    alabelstyle: LabelSpec;
    /** label for b axis */
    blabel: string;
    /** label style for a axis */
    blabelstyle: LabelSpec;
    /** label for c axis */
    clabel: string;
    /** label style for a axis */
    clabelstyle: LabelSpec;
}

export interface ArrowSpec {
    start: Vector3;
    end: Vector3;
    radius: number;
    color: ColorSpec;
    hidden: boolean;
    /** ratio of arrow base to cylinder (1.618034 default) */
    radiusRatio: number;
    /** relative position of arrow base (0.618034 default) */
    mid: number;
    /** position of arrow base in length units, if negative positioned from
     * end instead of start. Overrides mid. */
    midpos: number;
}

export interface LabelSpec {
    /** font name, default sans-serif */
    font: string;
    /** height of text, default 18 */
    fontSize: number;
    /** font color, default white */
    fontColor: ColorSpec;
    /** font opacity, default 1 */
    fontOpacity: number;
    /** line width of border around label, default 0 */
    borderThickness: number;
    /** color of border, default backgroundColor */
    borderColor: ColorSpec;
    /** color of border */
    borderOpacity: string;
    /** color of background, default black */
    backgroundColor: ColorSpec;
    /** opacity of background, default 1 */
    backgroundOpacity: string;
    /** x, y, z coordinates for label */
    position: Vector3;
    /** always put labels in front of model */
    inFront: boolean;
    /** show background rounded rectangle, default true */
    showBackground: boolean;
    /** sets the label to change with the model when zooming */
    fixed: boolean;
    /** An element to draw into the label. */
    backgroundImage: CanvasImageSource;
    /** how to orient the label w/respect to position: topLeft (default),
     * topCenter, topRight, centerLeft, center, centerRight, bottomLeft,
     * bottomCenter, bottomRight */
    alignment: string;
    /** if set, only display in this frame of an animation */
    frame: number;
}

export interface ShapeSpec {
    /** solid color */
    color: ColorSpec;
    /** transparency */
    alpha: number;
    /** draw as wireframe, not surface */
    wireframe: boolean;
    /** if true, do not display object */
    hidden: boolean;
    /** width of line for wireframe rendering. No longer supported by most browsers */
    linewidth: number;
    /** if true, user can click on object to trigger callback */
    clickable: boolean;
    /** function to call on click */
    callback: () => void;
    /** if set, only display in this frame of an animation */
    frame: number;
}

export type ClickSphereStyleSpec = unknown;

export declare class GLViewer {
    constructor(element: HTMLElement, config?: Partial<ViewerSpec>);

    public addModel(data?: string, format?: string, options?: Partial<ParserOptionsSpec>): GLModel;

    public render(): void;
    public removeModel(model: GLModel): void;

    public getFrame(): number;

    public setStyle(sel: Partial<AtomSelectionSpec>, style: Partial<AtomStyleSpec>): void;
    public selectedAtoms(sel: Partial<AtomSelectionSpec>): Array<Partial<AtomSpec>>;

    public createModelFrom(sel: Partial<AtomSelectionSpec>, extract?: boolean): GLModel;

    public removeAllModels(): void;

    public replicateUnitCell(A: number, B?: number, C?: number, model?: GLModel): void;

    public setClickable(
        sel: Partial<AtomSelectionSpec>,
        clickable: boolean,
        callback: (
            atom: Partial<AtomSpec>,
            viewer: GLViewer,
            event: unknown,
            container: unknown
        ) => void
    ): void;

    public zoomTo(
        sel?: Partial<AtomSelectionSpec>,
        animationDuration?: number,
        fixedPath?: boolean
    ): void;

    public zoom(factor?: number, animationDuration?: number, fixedPath?: boolean): void;

    public setSlab(near: number, far: number): void;
    public enableFog(enable: boolean): void;
    public resize(): void;

    public addUnitCell(model?: GLModel, spec?: Partial<UnitCellStyleSpec>): void;
    public removeUnitCell(model?: GLModel): void;

    public spin(axis?: 'x' | 'y' | 'z' | 'vx' | 'vy' | 'vz' | false): void;

    public getView(): View;
    public setView(view: View): void;
    public setCameraParameters(parameters: Partial<{ fov: number; z: number }>): void;

    public addArrow(spec: Partial<ArrowSpec>): GLShape;
    public removeShape(shape: GLShape): void;

    public addLabel(
        text: string,
        options?: Partial<LabelSpec>,
        sel?: Partial<AtomSelectionSpec>,
        noshow?: boolean
    ): Label;
    public removeLabel(label: Label): void;

    public pngURI(): string;

    public setHoverable(
        sel: Partial<AtomSelectionSpec>,
        hoverable: bool,
        hover_callback: (
            atom: AtomSpec,
            viewer: GLViewer,
            event: unknown,
            container: HTMLElement
        ) => void,
        unhover_callback: (atom: AtomSpec) => void
    ): void;

    public setHoverDuration(duration?: number): void;
}

/** [ pos.x, pos.y, pos.z, rotationGroup.position.z, q.x, q.y, q.z, q.w ] */
export type View = [number, number, number, number, number, number, number, number];

export declare class GLModel {
    public setStyle(
        sel: Partial<AtomSelectionSpec>,
        style: Partial<AtomStyleSpec>,
        add?: boolean
    ): void;

    public addAtoms(newatoms: Array<Partial<AtomSpec>>): void;

    public setCrystData(
        a: number,
        b: number,
        c: number,
        alpha: number,
        beta: number,
        gamma: number
    ): void;

    public setCrystMatrix(matrix: Matrix3): void;

    public selectedAtoms(sel: Partial<AtomSelectionSpec>): Array<Partial<AtomSpec>>;
    public removeAtoms(atoms: Array<Partial<AtomSpec>>): Array<Partial<AtomSpec>>;

    public addAtomSpecs(customAtomSpecs: string[]): void;

    public getInternalState(): { atoms: AtomSpec[]; frames: unknown[] };
    public setInternalState(state: { atoms: AtomSpec[]; frames: unknown[] }): void;
}

export declare class GLShape {}
export declare class Label {}

export declare function createViewer(
    element: HTMLElement | string,
    config?: Partial<ViewerSpec>,
    shared_viewer_resources?: unknown
): GLViewer | undefined;

interface PresetElementColors {
    defaultColor: ElementColors;
    Jmol: ElementColors;
    rasmol: ElementColors;
    defaultColors: ElementColors;
    greenCarbon: ElementColors;
    cyanCarbon: ElementColors;
    magentaCarbon: ElementColors;
    yellowCarbon: ElementColors;
    whiteCarbon: ElementColors;
    orangeCarbon: ElementColors;
    purpleCarbon: ElementColors;
    blueCarbon: ElementColors;
}

export declare const elementColors: PresetElementColors;

interface PresetSpriteAlignment {
    topLeft: Vector2;
    topCenter: Vector2;
    topRight: Vector2;
    centerLeft: Vector2;
    center: Vector2;
    centerRight: Vector2;
    bottomLeft: Vector2;
    bottomCenter: Vector2;
    bottomRight: Vector2;
}

export declare const SpriteAlignment: PresetSpriteAlignment;

export declare class Vector3 {
    constructor(x: number, y: number, z: number);
    public normalize(): Vector3;
}

export declare class Vector2 {
    constructor(x: number, y: number);
}

export declare class Matrix3 {
    // prettier-ignore
    constructor(
        ax: number, ay: number, az: number,
        bx: number, by: number, bz: number,
        cx: number, cy: number, cz: number
    );
}

export declare class Sphere {
    constructor(center: Vector3, radius: number);
}
