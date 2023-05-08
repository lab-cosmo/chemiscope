// import { StructureOptions } from '../../src/structure/options';

// import { assert } from 'chai';

// const DUMMY_CALLBACK = () => {
//     return { top: 0, left: 0 };
// };

// function createShadowChild(): HTMLElement {
//     const shadowElement = document.createElement('div');
//     const shadow = shadowElement.attachShadow({ mode: 'open' });

//     const element = document.createElement('div');
//     shadow.appendChild(element);

//     return element;
// }

// let KARMA_INSERTED_HTML: string;

// describe('StructureOptions', () => {
//     before(() => {
//         // store karma's default HTML
//         KARMA_INSERTED_HTML = document.body.innerHTML;
//     });

//     it('can remove itself from DOM', () => {
//         const root = createShadowChild();

//         const options = new StructureOptions(root, DUMMY_CALLBACK);
//         assert(root.innerHTML !== '');
//         assert(document.body.innerHTML !== KARMA_INSERTED_HTML);

//         options.remove();
//         assert(document.body.innerHTML === KARMA_INSERTED_HTML);
//         assert(root.innerHTML === '');
//     });
// });
