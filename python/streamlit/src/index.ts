// streamlit-component/src/index.ts
import {
  Streamlit,
  RenderData,
} from "streamlit-component-lib";

const ROOT_ID = "chemiscope-root";

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
  root.innerHTML = `<div class="container-fluid">
      <div class="row">
          <div class="col-md-6" style="padding: 0">
              <div class="ratio ratio-1x1">
                  <div id="chemiscope-meta" style="z-index: 10"></div>
                  <div id="chemiscope-map" style="position: absolute"></div>
              </div>
          </div>

          <div class="col-md-6" style="padding: 0">
              <div class="ratio ratio-5x7">
                  <div>
                      <div id="chemiscope-structure" style="height:450px"></div>
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

