function showPopup(type, data) {
    // Create new popup container dynamically
    const popup = document.createElement('div');
    popup.classList.add('popup', `popup-${type}`);

    popup.innerHTML = `
        <div class="popup-header">
            <span id="popup-title">Title</span>
            <span class="popup-close" onclick="this.closest('.popup').remove()">×</span>
        </div>
        <div class="popup-content" id="popup-content">
        Content
        </div>
    `;

    document.body.appendChild(popup);

    const title = popup.querySelector('#popup-title');
    const content = popup.querySelector('#popup-content');

    if (type === 'node') {
        const cid = data.id;
        const compound = cidToCompound.get(cid);
        const compoundProperties = cidToCompoundProperties.get(cid);
        title.textContent = compound ? compound.cmpdname : 'Node Info';
        
        let html = '<div class="compound-structure">';
        html += `<img src="data/structures/${cid}.svg" alt="${compound.cmpdname}" onerror="this.style.display='none'">`;
        html += '</div>';

        console.log(cidToHazards.get(cid));
        if (cidToHazards.has(cid)) {
            console.log("lul")
            const compoundHazards = cidToHazards.get(cid);
            html += '<div class="compound-hazards">'
            if (compoundHazards.nfpa) {

                let createNfpaDiamond = (health, flammability, instability) => {
                    return `
                    <div class="nfpa-diamond">
                    <div class="nfpa-diamond-grid">

                        <div class="nfpa-diamond-item nfpa-flammability">
                        <span class="nfpa-diamond-item-value">${flammability}</span>
                        </div>

                        <div class="nfpa-diamond-item nfpa-stability">
                        <span class="nfpa-diamond-item-value">${instability}</span>
                        </div>

                        <div class="nfpa-diamond-item nfpa-health">
                        <span class="nfpa-diamond-item-value">${health}</span>
                        </div>

                        <div class="nfpa-diamond-item nfpa-special">
                        <span class="nfpa-diamond-item-value"></span>
                        </div>

                    </div>
                </div>`
                };


                const nfpa = compoundHazards.nfpa;
                const health = 'health' in nfpa ? nfpa.health.toString() : "";
                const flammability = 'flammability' in nfpa ? nfpa.flammability.toString() : "";
                const instability = 'instability' in nfpa ? nfpa.instability.toString() : "";
                
                html += createNfpaDiamond(health, flammability, instability);
            }

            if (compoundHazards.pictograms) {
                html += '<div class="ghs-pictograms">';

                pictograms = compoundHazards.pictograms;
                for (ghs of pictograms) {
                    html += `<img src="data/assets/ghs_pictograms/${ghs}.svg">`
                }
                
                html += '</div>';
            }

            html += '</div>';
        }

        html += '<div class="compound-info">';
        for (const entry of compoundProperties) {
            const label = entry.property;
            const value = entry.value.split('\n').join('<br>');
            html += '<div class="compound-info-row">';
            html += `<span class="compound-info-label table-info">${label}</span>`;
            html += '<span class="compound-info-value table-info">'
            if (label.toLowerCase().includes('pubchem') && label.toLowerCase().includes('cid'))
                html += `<a href="https://pubchem.ncbi.nlm.nih.gov/compound/${value}" target="_blank">${value}</a>`;
            else if (label.toLowerCase().includes('wikipedia'))
                html += `<a href="${value}" target="_blank">${value.split('/').at(-1).replace('_', ' ')}</a>`;
            else
                html += value;
            html += '</span>';
            html += '</div>';
        }
        html += '</div>';
        
        const description = cidToDescription.get(cid);
        if (description) {
            html += '<div class="compound-info">';
            html += '<span class="compound-info-label description">Description:</span>';
            const description_html = '<p>' + description.split('\n\n').join('</p><p>') + '</p>';
            html += `<span class="compound-info-value description">${description_html}</span>`;
            html += '</div>';
            html += '</div>';

        }
        
        content.innerHTML = html;
    } else if (type === 'edge') {
        const sourceId = typeof data.source === 'object' ? data.source.id : data.source;
        const targetId = typeof data.target === 'object' ? data.target.id : data.target;
        const sourceComp = cidToCompound.get(sourceId);
        const targetComp = cidToCompound.get(targetId);
        const sourceName = sourceComp ? sourceComp.cmpdname : sourceId;
        const targetName = targetComp ? targetComp.cmpdname : targetId;
        
        title.textContent = 'Reactions';
        
        const edgeStr = `${sourceId}-${targetId}`;
        const reverseEdgeStr = `${targetId}-${sourceId}`;
        const forwardReactionIDs = edgeToReactionID.get(edgeStr) || [];
        const reverseReactionIDs = data.type === 'bi' ? (edgeToReactionID.get(reverseEdgeStr) || []) : [];
        
        let html = '';

        const generateReactionParticipants = (participants, balanced) => {
            let reactants_htmls = [];
            for (const entry of participants) {
                const cid = entry.cid;
                const compound = cidToCompound.get(cid);
                const name = compound.cmpdname;
                let curr_html = ""
                curr_html += `<div class="reaction-participant">`;
                curr_html += `<img class="reaction-svg-image" src="data/structures/${cid}.svg" alt="${name}" data-cid="${cid}">`;
                curr_html += `<span class="reaction-participant-name" title="${name}">${name}</span>`;
                curr_html += '</div>';
                
                if (balanced) {
                    const coeff = entry.coeff;
                    curr_html = `<div class="reaction-participant-balanced"><span class="reaction-participant-coeff">${coeff}</span>${curr_html}</div>`
                }
                reactants_htmls.push(curr_html)
            }
            return reactants_htmls.join('<span class="reaction-sep reaction-sep-plus">+</span>')
        };

        const generateReactionItem = (rid) => {
            const getConfidenceClass = (confidence) => {
                const classPrefix = 'reaction-confidence-'
                if (confidence < 0.5)
                    return classPrefix + 'low';
                else if (confidence < 0.7)
                    return classPrefix + 'medium';
                else
                    return classPrefix + 'high';
            };

            item_html = "";
            const reaction = RIDToReaction.get(rid);
            if (reaction) {
                const balanced = reaction.balanced;
                console.log(balanced)
                item_html += '<div class="reaction-item">';
                item_html += `<div class="reaction-equation">`;
                item_html += generateReactionParticipants(reaction.reagents, balanced);
                item_html += '<span class="reaction-sep reaction-sep-arrow">→</span>';
                item_html += generateReactionParticipants(reaction.products, balanced);
                item_html += '</div>';

                const description = ridToDescription.get(rid);
                if (description) {
                    item_html += '<div class="reaction-description">';
                    item_html += '<strong>Description:</strong> ';
                    item_html += description;
                    item_html += '</div>';
                }

                const confidence = reaction.confidence;
                if (confidence) {
                    const confidenceClass = getConfidenceClass(confidence);
                    item_html += `<span class="reaction-confidence ${confidenceClass}" title="This reaction is generated automatically and this is its validation score.">${confidence.toFixed(2)}</span>`;
                } else if (reaction.source == 'ord') {
                    item_html += `<img src="data/assets/ord.svg" class="reaction-ord-label" title="Sourced from Open Reaction Database">`
                }

                item_html += '</div>';
            }

            return item_html;
        };

        const generateReactionList = (rids) => {
            let reactionsHtml = "";
            if (rids.length > 0) {
                rids.forEach(rid => {
                    reactionsHtml += generateReactionItem(rid);
                });
            } else {
                reactionsHtml += '<div class="reaction-item">No reactions found</div>';
            }

            return reactionsHtml;
        };
        
        if (data.type === 'bi') {
            html += '<div class="popup-split">';
            
            html += '<div class="popup-direction">';
            html += `<h4> ${sourceName} → ${targetName}</h4>`;
            html += generateReactionList(forwardReactionIDs);
            html += '</div>';
            
            html += '<div class="popup-direction">';
            html += `<h4> ${targetName} → ${sourceName}</h4>`;
            html += generateReactionList(reverseReactionIDs);
            html += '</div>';
            
            html += '</div>';
        } else {
            html += '<div class="popup-direction">';
            html += `<h4> ${sourceName} → ${targetName}</h4>`;
            html += generateReactionList(forwardReactionIDs);
            html += '</div>';
        }
        
        content.innerHTML = html;
    }

    popup.style.display = 'block';
    popup.style.position = 'absolute';
    popup.style.visibility = 'visible';
    popup.style.opacity = '1';

    const rect = popup.getBoundingClientRect();
    popup.style.left = (window.innerWidth / 2 - rect.width / 2) + 'px';
    popup.style.top = (window.innerHeight / 2 - rect.height / 2) + 'px';

    makeDraggable(popup);

    document.querySelectorAll('.reaction-svg-image').forEach(el => {
        el.onclick = () => {
            const cid = Number(el.dataset.cid);
            showCompoundInfoPopup(cid);
        };
    });
}


function showCompoundInfoPopup(cid) {
    const name = cidToCompound.get(cid).cmpdname;
    showPopup('node', {id: cid, name: name});
}