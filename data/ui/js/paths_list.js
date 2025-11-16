
class PathsList {
    constructor(compounds) {
        this.selectedLevel = 1;
        this.sourceCids = new Set();
        this.targetCids = new Set();
        this.compounds = compounds
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
        const cidPrefix = document.createTextNode('CID: ');
        itemExtra.appendChild(cidPrefix);
        itemExtra.appendChild(cidLink);

        if (compound.wiki) {
            const wikiLink = document.createElement('a');
            wikiLink.href = compound.wiki;
            wikiLink.textContent = `${compound.wiki.split('/').at(-1).replace('_', ' ')}`;
            wikiLink.target = '_blank';
            
            const sep_wiki_prefix = document.createTextNode(' | Wiki: ');
            itemExtra.appendChild(sep_wiki_prefix);
            itemExtra.appendChild(wikiLink);
        }
        
        itemInfo.appendChild(itemName);
        itemInfo.appendChild(itemExtra);

        return itemInfo;
    }


    selectLevel(level) {
        this.selectedLevel = level;
        document.getElementById('level1').classList.toggle('selected', level === 1);
        document.getElementById('level2').classList.toggle('selected', level === 2);
    }


    #createPathItem(compound) {
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

        removeBtn.addEventListener('click', () => {
            this.removeCidAndRenderLists(cid);
        });

        item.appendChild(iconContainer);
        item.appendChild(itemInfo);
        item.appendChild(removeBtn);

        item.addEventListener('contextmenu', (e) => {
            e.preventDefault();
            this.removeCidAndRenderLists(cid);
        });

        return item;
    }


    renderPathList() {
        const level1Items = document.querySelector('#level1 .level-items');
        const level2Items = document.querySelector('#level2 .level-items');
        level1Items.innerHTML = '';
        level2Items.innerHTML = '';

        this.sourceCids.forEach(cid => {
            const comp = this.compounds.get(cid);
            if (comp) {
                const item = this.#createPathItem(comp);
                level1Items.appendChild(item);
            }
        });

        this.targetCids.forEach(cid => {
            const comp = this.compounds.get(cid);
            if (comp) {
                const item = this.#createPathItem(comp);
                level2Items.appendChild(item);
            }
        });
    }


    removeCidAndRenderLists(cid) {
        if (this.selectedLevel === 1) {
            this.sourceCids.delete(cid);
        } else {
            this.targetCids.delete(cid);
        }
        this.renderPathList();
    }


    addCidAndRenderLists(cid) {
        if (this.selectedLevel === 1) {
            this.sourceCids.add(cid);
        } else {
            this.targetCids.add(cid);
        }
        this.renderPathList();
    }
}
