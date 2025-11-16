// streamlit-component/src/index.ts
import {
  Streamlit,
  RenderData,
} from "streamlit-component-lib";

const ROOT_ID = "chemiscope-root";
let visualizer: any | null = null;
let indexer: any | null = null;
let visualizerLoaded = false;
let lastSelection: number | null = null;
let lastReportedSelection: number | null = null; // CS -> ST
let originalMapOnselect: any | null = null;
let originalStructOnselect: any | null = null;
let originalInfoOnchange: any | null = null;

// Get the Chemiscope object from window
function getChemiscope(): any | null {
  const cs = (window as any).Chemiscope;
  if (!cs) {
    console.error(
      "window.chemiscope is not defined. " +
        "Did chemiscope.min.js load in index.html?"
    );
    return null;
  }
  if (!cs.DefaultVisualizer) {
    console.error("window.chemiscope.DefaultVisualizer is missing:", cs);
    return null;
  }
  return cs;
}

function getOrCreateRoot(): HTMLDivElement {
  console.log("getting or creating root");
  let root = document.getElementById(ROOT_ID) as HTMLDivElement | null;
  if (!root) {
    root = document.createElement("div");
    root.id = ROOT_ID;
    root.style.width = "100%";
    document.body.appendChild(root);
  }
  return root;
}

// Apply a selection coming from Streamlit (structure index)
function applySelectionFromStructure(structureIndex: number): void {
  if (!visualizer || !indexer) return;

  // We interpret the Streamlit slider as controlling per-structure selection.
  // Use the indexer in "structure" target mode to get full Indexes.
  const indexes = indexer.fromStructureAtom("structure", structureIndex);
  if (!indexes) {
    console.warn("No environment for structure index", structureIndex);
    return;
  }

  try {
    console.log("selecting structure from ST", indexes)
    // 1) Move the active marker on the map
    visualizer.map.select(indexes);

    // 2) Mirror what a user click does: notify the rest of the system
    //    (DefaultVisualizer wired map.onselect to update info + structure)
    console.log("maponselect", originalMapOnselect)
    if (originalMapOnselect) {
      originalMapOnselect(indexes); // uses the unwrapped version to avoid infinite loops
    }
  } catch (err) {
    console.error("Error applying selection to chemiscope from slider:", err);
  }
}

function reportSelectionToStreamlit(structureIndex: number): void {
  console.log("reporting selection to Streamlit:", structureIndex); 
  if (structureIndex === lastReportedSelection) {
    return; // nothing new
  }
  lastReportedSelection = structureIndex;

  try {
    Streamlit.setComponentValue(structureIndex);
  } catch (err) {
    console.error("Error sending selection back to Streamlit:", err);
  }
}

function installReverseSyncCallbacks(): void {
  if (!visualizer) return;
  console.log("installing callbacks")

  const sendFromIndexes = (indexes: any) => {
    if (!indexes || typeof indexes.structure !== "number") return;
    const s = indexes.structure as number;
    // update our forward-selection cache as well
    lastSelection = s;
    reportSelectionToStreamlit(s);
  };

  // Wrap map.onselect
  originalMapOnselect = visualizer.map.onselect;
  visualizer.map.onselect = (indexes: any) => {
    console.log("map.onselect called with indexes:", indexes);
    if (typeof originalMapOnselect === "function") {
      originalMapOnselect(indexes);
    }
    sendFromIndexes(indexes);
  };

  // Wrap structure.onselect
  originalStructOnselect = visualizer.structure.onselect;
  visualizer.structure.onselect = (indexes: any) => {
    console.log("structure.onselect called with indexes:", indexes);
    if (typeof originalStructOnselect === "function") {
      originalStructOnselect(indexes);
    }
    sendFromIndexes(indexes);
  };

  // Wrap info.onchange
  originalInfoOnchange = visualizer.info.onchange;
  visualizer.info.onchange = (indexes: any) => {
    console.log("info.onchange called with indexes:", indexes);
    if (typeof originalInfoOnchange === "function") {
      originalInfoOnchange(indexes);
    }
    sendFromIndexes(indexes);
  };
}



function onRender(event: Event): void {
  const data = (event as CustomEvent<RenderData>).detail;
  const args = data.args as {
    dataset: any;
    height?: number;
    selected_index?: number;
  };

  const dataset = args.dataset;
  const height = typeof args.height === "number" ? args.height : 600;
  const selectedIndex = typeof args.selected_index === "number"
    ? args.selected_index
    : null;

  if (!dataset) {
    return;
  }

  const Chemiscope = getChemiscope();
  if (!Chemiscope) {
    return;
  }

  // First time: load the visualizer
  if (!visualizerLoaded) {
    visualizerLoaded = true;

    const root = getOrCreateRoot();
      root.innerHTML = `<div class="container-fluid" style="padding: 0; background-color: white; opacity: 1;">
          <div class="row">
              <div class="col-md-6" style="padding: 0">
                  <div class="ratio ratio-1x1">
                      <div id="chemiscope-meta" style="z-index: 10; height: 50px"></div>
                      <div id="chemiscope-map" style="position: absolute"></div>
                  </div>
              </div>

              <div class="col-md-6" style="padding: 0">
                  <div class="ratio ratio-5x7">
                      <div>
                          <div id="chemiscope-structure" style="height:420px"></div>
                          <div id="chemiscope-info"></div>
                      </div>
                  </div>
              </div>
          </div>
      </div>`;

      const config = {
          map: 'chemiscope-map',
          info: 'chemiscope-info',
          meta: 'chemiscope-meta',
          structure: 'chemiscope-structure',
      };
          

      Chemiscope.DefaultVisualizer.load(config as any, dataset as any)
        .then((v: any) => {
          visualizer = v;

          // Build our own EnvironmentIndexer from the current dataset
          const ds = v.dataset(true); // get structures resolved
          indexer = new Chemiscope.EnvironmentIndexer(
            ds.structures,
            ds.environments
          );

          installReverseSyncCallbacks(); 

          // If we already have a pending selection, apply it
          if (selectedIndex !== null) {
            lastSelection = selectedIndex;
            applySelectionFromStructure(selectedIndex);
          }
        })
        .catch((err: unknown) => {
          console.error("Error loading chemiscope visualizer:", err);
        });
  } else {
    // Subsequent renders: just apply the new selection, if any
    if (
      selectedIndex !== null &&
      selectedIndex !== lastSelection &&
      visualizer
    ) {
      lastSelection = selectedIndex;
      applySelectionFromStructure(selectedIndex);
    }
  }

  Streamlit.setFrameHeight(height + 40);
}

Streamlit.events.addEventListener(Streamlit.RENDER_EVENT, onRender);
Streamlit.setComponentReady();


/*
function onRender(event: Event): void {
  console.log("rendering")
  const data = (event as CustomEvent<RenderData>).detail;
  const args = data.args as { dataset: any; height?: number };

  const dataset = args.dataset;
  const height = typeof args.height === "number" ? args.height : 600;

  if (!dataset) {
    return; // nothing to render yet
  }

  const Chemiscope = getChemiscope();
  if (!Chemiscope) {
    return; // error already logged
  }

  const root = getOrCreateRoot();
  root.innerHTML = `<div class="container-fluid" style="padding: 0; background-color: white; opacity: 1;">
      <div class="row">
          <div class="col-md-6" style="padding: 0">
              <div class="ratio ratio-1x1">
                  <div id="chemiscope-meta" style="z-index: 10; height: 50px"></div>
                  <div id="chemiscope-map" style="position: absolute"></div>
              </div>
          </div>

          <div class="col-md-6" style="padding: 0">
              <div class="ratio ratio-5x7">
                  <div>
                      <div id="chemiscope-structure" style="height:420px"></div>
                      <div id="chemiscope-info"></div>
                  </div>
              </div>
          </div>
      </div>
  </div>`;

  const config = {
      map: 'chemiscope-map',
      info: 'chemiscope-info',
      meta: 'chemiscope-meta',
      structure: 'chemiscope-structure',
  };

  console.log("loading chemiscope visualizer");

  // Load chemiscope visualizer
  Chemiscope.DefaultVisualizer.load(config as any, dataset as any).catch(
    (err: unknown) => {
      // eslint-disable-next-line no-console
      console.error("Error loading chemiscope visualizer:", err);
    },
  );

  // Adjust iframe height
  Streamlit.setFrameHeight(height + 40);
}

// register render handler
Streamlit.events.addEventListener(Streamlit.RENDER_EVENT, onRender);

// let Streamlit know we're ready
Streamlit.setComponentReady();
*/
