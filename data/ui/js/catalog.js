class Catalog extends PathsList {
    constructor(compounds, popup) {
        super(compounds);
        this.popup = popup;
        this.results = [];
        this.isEnd = false;
        this.query = null;
        this.sortingOrder = "none";
        this.sortingDirection = "asc"
        this.currentPage = 0;
        this.searchApiEndpoint = "/api/search";
        this.loaded = null;

        this.loadingHtml = document.getElementById('results-container').innerHTML;
    }

    [Symbol.iterator]() {
        let index = 0;
        const results = this.results;

        return {
            next() {
                if (index < results.length) {
                    return { value: results[index++], done: false };
                } else {
                    return { done: true };
                }
            }
        };
    }

    setQuery(query = "") {
        if(typeof query !== "string")
            throw new Error("Query must be string");

        this.query = query.trim();
        this.results = [];
        this.isEnd = false;
        this.currentPage = 0;

        this.loadMoreAndDisplay();
    }

    changeSortingOrder(sortingOrder, sortingDirection) {
        this.sortingOrder = sortingOrder;
        this.sortingDirection = sortingDirection;
        this.setQuery(this.query);
    }

    loadMoreAndDisplay() {
        this.loaded = this.#loadMore();
        this.displayResults();
    }

    async #loadMore() {
        if (this.isEnd) return;

        const url = new URL(this.searchApiEndpoint, window.location.origin);
        const params = url.searchParams;

        params.append("sorting_order", this.sortingOrder);
        params.append("sorting_direction", this.sortingDirection);
        params.append("page", this.currentPage.toString());

        if (this.query !== "") {
            params.append("q", this.query);
        }

        try {
            const response = await fetch(url);
            if (!response.ok) throw new Error(`HTTP error during compounds fetch! status: ${response.status}`);

            const data = await response.json();
            this.results.push(...data.cids);
            this.isEnd = data.is_end;
            this.currentPage += 1;

            return this.compounds.loadCompounds(data.cids);

        } catch (error) {
            console.error('Error loading results:', error);
            return false;
        }
    }

    #createResultItem(cid) {
        const compound = this.compounds.get(cid);
        const item = document.createElement('div');
        item.className = 'result-item';
        
        const checkbox = document.createElement('input');
        checkbox.type = 'checkbox';
        checkbox.className = 'item-check';
        checkbox.checked = selectedCIDs.has(compound.cid);
        checkbox.addEventListener('change', function() {
            if (this.checked) {
                selectedCIDs.add(compound.cid);
                // Updates graph adding edges attached to this node
                updateGraph(compound.cid);
            } else {
                // Removes node with attached edges
                removeNode(compound.cid);
                // Updates graph without updating edges
                updateGraph();
            }
        });
        
        const iconContainer = document.createElement('div');
        iconContainer.className = 'item-icon loading';
        iconContainer.textContent = 'Loading...';
        iconContainer.dataset.cid = `${compound.cid}`;
        iconContainer.addEventListener('click', (e) => {
            const target = e.currentTarget;
            const cid = Number(target.dataset.cid);
            showCompoundInfoPopup(cid);
        });
        
        loadStructureSVG(cid, iconContainer);
        
        const itemInfo = this.createItemInfo(compound);

        const addItemToPathList = () => {
            this.addCidAndRenderLists(compound.cid);
        };

        const addBtn = document.createElement('button');
        addBtn.className = 'compound-item-btn add-item';
        addBtn.textContent = '+';
        addBtn.addEventListener('click', () => {
            addItemToPathList();
        });
        
        item.appendChild(checkbox);
        item.appendChild(iconContainer);
        item.appendChild(itemInfo);
        item.appendChild(addBtn);

        item.addEventListener('contextmenu', (e) => {
            e.preventDefault();
            addItemToPathList();
        });
        
        return item;
    }

    #addLoadMoreButton() {
        const resultsContainer = document.getElementById('results-container');
        const loadMoreBtn = document.createElement('button');
        loadMoreBtn.className = 'load-more-btn';
        loadMoreBtn.textContent = 'Load More Results';
        loadMoreBtn.onclick = this.loadMoreAndDisplay.bind(this);
        resultsContainer.appendChild(loadMoreBtn);
    }

    async displayResults() {
        const resultsContainer = document.getElementById('results-container');
        resultsContainer.innerHTML = this.loadingHtml;

        if(this.loaded === null || !await this.loaded) {
            resultsContainer.innerHTML = '<div class="no-results">Error occured during compounds loading</div>';
            return;
        }
        if (this.results.length === 0) {
            resultsContainer.innerHTML = '<div class="no-results">No compounds found</div>';
            return;
        }

        const existingBtn = resultsContainer.querySelector('.load-more-btn');
        if (existingBtn) {
            existingBtn.remove();
        }
        
        resultsContainer.innerHTML = '';
        for(const cid of this.results) {
            const resultItem = this.#createResultItem(cid);
            resultsContainer.appendChild(resultItem);
        }

        if (!this.isEnd) {
            this.#addLoadMoreButton();
        }
    }


    showPopupNode(cid) {
        const compound = compounds.get(cid);
        const name = compound.name;
        
        let html = '<div class="compound-structure">';
        html += `<img src="${compounds.getSvgPath(cid)}" alt="${name}" onerror="this.style.display='none'">`;
        html += '</div>';

        if (compound.nfpa || compound.pictograms) {
            html += '<div class="compound-hazards">'
            if (compound.nfpa) {

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


                const nfpa = compound.nfpa;
                const health = 'health' in nfpa ? nfpa.health : "";
                const flammability = 'flammability' in nfpa ? nfpa.flammability : "";
                const instability = 'instability' in nfpa ? nfpa.instability : "";
                
                html += createNfpaDiamond(health, flammability, instability);
            }

            if (compound.pictograms) {
                html += '<div class="ghs-pictograms">';

                const pictograms = compound.pictograms;
                for (const ghs of pictograms) {
                    html += `<img src="/assets/ghs_pictograms/${ghs}.svg">`
                }
                
                html += '</div>';
            }

            html += '</div>';
        }

        html += '<div class="compound-info compound-info-general">';
        for (const entry of compound.properties) {
            const label = entry.property_name;
            const value = entry.property_value.split('\n').join('<br>');
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
        console.log(compound);
        const description = compound.description;
        if (description) {
            html += '<div class="compound-info compound-info-description">';
            const description_html = '<p>' + description.split('\n\n').join('</p><p>') + '</p>';
            html += `<span class="compound-info-value description">${description_html}</span>`;
            html += '</div>';
            html += '</div>';
        }

        this.popup.showPopup(html, name, 'popup-node');
    }
}