import { StructureOptions } from '../../src/structure/options';
import { GUID } from '../../src/utils';

import { assert } from 'chai';

const DUMMY_CALLBACK = () => {
    return { top: 0, left: 0 };
};

function createShadowChild(): HTMLElement {
    const shadowElement = document.createElement('div');
    const shadow = shadowElement.attachShadow({ mode: 'open' });

    const element = document.createElement('div');
    shadow.appendChild(element);

    return element;
}

function traverseDOM(element: Element, callback: (element: Element) => void) {
    if (element.children) {
        const elements = element.children;
        for (let i = 0; i < elements.length; i++) {
            traverseDOM(elements[i], callback);
        }
    }
    callback(element);
}

let KARMA_INSERTED_HTML: string;

describe('StructureOptions', () => {
    before(() => {
        // store karma's default HTML
        KARMA_INSERTED_HTML = document.body.innerHTML;
    });

    it('can remove itself from DOM', () => {
        const root = createShadowChild();

        const options = new StructureOptions(root, 'this-is-my-id' as GUID, DUMMY_CALLBACK);
        assert(root.innerHTML !== '');
        assert(document.body.innerHTML !== KARMA_INSERTED_HTML);

        options.remove();
        assert(document.body.innerHTML === KARMA_INSERTED_HTML);
        assert(root.innerHTML === '');
    });

    it('has a unique id in the page', () => {
        const root = createShadowChild();

        const guid = 'this-is-my-id' as GUID;
        const options = new StructureOptions(root, guid, DUMMY_CALLBACK);
        traverseDOM(document.body, (element) => {
            if (element.id) {
                assert(element.id.includes(guid));
            }
        });

        options.remove();
    });
});
