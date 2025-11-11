class ItemsList {

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
}