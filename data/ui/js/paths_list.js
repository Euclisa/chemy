class PathsList {
    constructor(compounds) {
        this.selectedLevel = 1;
        this.sourceCids = new Set();
        this.targetCids = new Set();
        this.compounds = compounds;
    }


    createItemInfo(compound) {
        const itemInfo = document.createElement('div');
        itemInfo.className = 'item-info';

        const itemName = document.createElement('div');
        itemName.className = 'item-name';
        itemName.textContent = compound.name || 'Unknown compound';

        const itemExtra = document.createElement('div');
        itemExtra.className = 'item-extra';

        const cidLink = document.createElement('a');
        cidLink.href = `https://pubchem.ncbi.nlm.nih.gov/compound/${compound.cid}`;
        cidLink.textContent = `${compound.cid}`;
        cidLink.target = '_blank';
        cidLink.rel = 'noopener noreferrer';
        itemExtra.appendChild(document.createTextNode('CID: '));
        itemExtra.appendChild(cidLink);

        if (compound.wiki) {
            const wikiLink = document.createElement('a');
            wikiLink.href = compound.wiki;
            wikiLink.textContent = `${compound.wiki.split('/').at(-1).replace(/_/g, ' ')}`;
            wikiLink.target = '_blank';
            wikiLink.rel = 'noopener noreferrer';

            itemExtra.appendChild(document.createTextNode(' | Wiki: '));
            itemExtra.appendChild(wikiLink);
        }

        itemInfo.appendChild(itemName);
        itemInfo.appendChild(itemExtra);

        return itemInfo;
    }


    hasCidInLevel(cid, level = this.selectedLevel) {
        return level === 1 ? this.sourceCids.has(cid) : this.targetCids.has(cid);
    }


    selectLevel(level) {
        this.selectedLevel = level;
        document.getElementById('level1').classList.toggle('selected', level === 1);
        document.getElementById('level2').classList.toggle('selected', level === 2);
    }


    #createPathItem(compound, level) {
        const cid = compound.cid;

        const item = document.createElement('div');
        item.className = 'result-item';
        item.dataset.cid = cid;

        const iconContainer = document.createElement('div');
        iconContainer.className = 'item-icon loading';
        iconContainer.textContent = 'Loading...';
        loadStructureSVG(cid, iconContainer);

        const itemInfo = this.createItemInfo(compound);

        const removeBtn = document.createElement('button');
        removeBtn.className = 'compound-item-btn remove-item';
        removeBtn.textContent = 'X';
        removeBtn.type = 'button';
        removeBtn.addEventListener('click', () => {
            this.removeCidAndRenderLists(cid, level);
        });

        item.appendChild(iconContainer);
        item.appendChild(itemInfo);
        item.appendChild(removeBtn);

        item.addEventListener('contextmenu', (e) => {
            e.preventDefault();
            this.removeCidAndRenderLists(cid, level);
        });

        return item;
    }


    renderPathList() {
        const level1Items = document.querySelector('#level1 .level-items');
        const level2Items = document.querySelector('#level2 .level-items');
        level1Items.innerHTML = '';
        level2Items.innerHTML = '';

        this.sourceCids.forEach((cid) => {
            const compound = this.compounds.get(cid);
            if (compound) {
                level1Items.appendChild(this.#createPathItem(compound, 1));
            }
        });

        this.targetCids.forEach((cid) => {
            const compound = this.compounds.get(cid);
            if (compound) {
                level2Items.appendChild(this.#createPathItem(compound, 2));
            }
        });
    }


    removeCidAndRenderLists(cid, level = this.selectedLevel) {
        const cids = level === 1 ? this.sourceCids : this.targetCids;
        cids.delete(cid);
        this.renderPathList();
    }


    addCidAndRenderLists(cid, level = this.selectedLevel) {
        const cids = level === 1 ? this.sourceCids : this.targetCids;
        cids.add(cid);
        this.renderPathList();
    }


    resetPathLists() {
        this.sourceCids.clear();
        this.targetCids.clear();
        this.selectLevel(1);
        this.renderPathList();
    }
}
