class Catalog extends PathsList {
    constructor(compounds, popup) {
        super(compounds);
        this.popup = popup;
        this.results = [];
        this.isEnd = false;
        this.query = "";
        this.sortingOrder = "none";
        this.sortingDirection = "asc";
        this.currentPage = 0;
        this.searchApiEndpoint = "/api/search";
        this.loaded = null;
        this.isLoading = false;
        this.loadingToken = 0;
        this.requestToken = 0;

        this.sortBtn = document.getElementById('catalog-sort-btn');
        this.sortBtnInitContent = this.sortBtn.innerHTML;
        document.querySelectorAll('.filter-menu li').forEach((item) => {
            item.addEventListener('click', () => {
                this.sortBtn.innerHTML = item.innerHTML;
                this.changeSortingOrder(item.dataset.sortOrder, item.dataset.sortDirection);
                this.refreshResults();
            });
        });

        this.searchButton = document.getElementById('search-button');
        this.searchInput = document.getElementById('search-input');
        this.handleSearch = () => {
            this.resetSortingOrder();
            this.setQuery(this.searchInput.value);
        };

        this.searchButton.addEventListener('click', this.handleSearch);
        this.searchInput.addEventListener('keydown', (e) => {
            if (e.key === 'Enter') {
                this.handleSearch();
            }
        });

        this.loadingHtml = document.getElementById('results-container').innerHTML;
    }


    setQuery(query = "") {
        if (typeof query !== "string") {
            throw new Error("Query must be string");
        }

        this.query = query.trim();
        this.results = [];
        this.isEnd = false;
        this.currentPage = 0;
        this.requestToken += 1;
        this.loadMoreAndDisplay(this.requestToken);
    }


    refreshResults() {
        this.setQuery(this.query);
    }


    changeSortingOrder(sortingOrder, sortingDirection) {
        this.sortingOrder = sortingOrder;
        this.sortingDirection = sortingDirection;
    }


    resetSortingOrder() {
        this.changeSortingOrder("none", "asc");
        this.sortBtn.innerHTML = this.sortBtnInitContent;
    }


    selectLevel(level) {
        super.selectLevel(level);
        this.syncSelectedLevelInputs();
    }


    addCidAndRenderLists(cid, level = this.selectedLevel) {
        super.addCidAndRenderLists(cid, level);
        this.syncSelectedLevelInputs();
    }


    removeCidAndRenderLists(cid, level = this.selectedLevel) {
        super.removeCidAndRenderLists(cid, level);
        this.syncSelectedLevelInputs();
    }


    resetPathLists() {
        super.resetPathLists();
        this.syncSelectedLevelInputs();
    }


    syncSelectedLevelInputs() {
        document.querySelectorAll('.item-check[data-cid]').forEach((input) => {
            const cid = Number(input.dataset.cid);
            input.checked = this.hasCidInLevel(cid);
        });
    }


    loadMoreAndDisplay(requestToken = this.requestToken) {
        if (this.isLoading && requestToken === this.loadingToken) {
            return;
        }

        this.isLoading = true;
        this.loadingToken = requestToken;
        this.loaded = this.#loadMore(requestToken);
        this.displayResults(requestToken);
    }


    async #loadMore(requestToken) {
        if (this.isEnd) {
            return { ok: true, stale: false };
        }

        const url = new URL(this.searchApiEndpoint, window.location.origin);
        url.searchParams.append("sorting_order", this.sortingOrder);
        url.searchParams.append("sorting_direction", this.sortingDirection);
        url.searchParams.append("page", this.currentPage.toString());

        if (this.query !== "") {
            url.searchParams.append("q", this.query);
        }

        try {
            const response = await fetch(url);
            if (!response.ok) {
                throw new Error(`HTTP error during compounds fetch! status: ${response.status}`);
            }

            const data = await response.json();
            if (requestToken !== this.requestToken) {
                return { ok: false, stale: true };
            }

            this.results.push(...data.cids);
            this.isEnd = data.is_end;
            this.currentPage += 1;

            const compoundsLoaded = await this.compounds.loadCompounds(data.cids);
            if (requestToken !== this.requestToken) {
                return { ok: false, stale: true };
            }

            return { ok: compoundsLoaded, stale: false };
        } catch (error) {
            console.error('Error loading results:', error);
            return { ok: false, stale: false };
        } finally {
            if (requestToken === this.loadingToken) {
                this.isLoading = false;
            }
        }
    }


    #createResultItem(cid) {
        const compound = this.compounds.get(cid);
        const item = document.createElement('div');
        item.className = 'result-item';

        const checkbox = document.createElement('input');
        checkbox.type = 'checkbox';
        checkbox.className = 'item-check';
        checkbox.dataset.cid = `${compound.cid}`;
        checkbox.checked = this.hasCidInLevel(compound.cid);
        checkbox.title = 'Add or remove this compound from the selected path list';
        checkbox.addEventListener('change', () => {
            if (checkbox.checked) {
                this.addCidAndRenderLists(compound.cid);
            } else {
                this.removeCidAndRenderLists(compound.cid);
            }
        });

        const iconContainer = document.createElement('div');
        iconContainer.className = 'item-icon loading';
        iconContainer.textContent = 'Loading...';
        iconContainer.dataset.cid = `${compound.cid}`;
        iconContainer.addEventListener('click', (e) => {
            const target = e.currentTarget;
            this.showPopupNode(Number(target.dataset.cid));
        });
        loadStructureSVG(cid, iconContainer);

        const itemInfo = this.createItemInfo(compound);
        const addBtn = document.createElement('button');
        addBtn.className = 'compound-item-btn add-item';
        addBtn.textContent = '+';
        addBtn.type = 'button';
        addBtn.addEventListener('click', () => {
            this.addCidAndRenderLists(compound.cid);
        });

        item.appendChild(checkbox);
        item.appendChild(iconContainer);
        item.appendChild(itemInfo);
        item.appendChild(addBtn);

        item.addEventListener('contextmenu', (e) => {
            e.preventDefault();
            this.addCidAndRenderLists(compound.cid);
        });

        return item;
    }


    #addLoadMoreButton() {
        const resultsContainer = document.getElementById('results-container');
        const loadMoreBtn = document.createElement('button');
        loadMoreBtn.className = 'load-more-btn';
        loadMoreBtn.type = 'button';
        loadMoreBtn.textContent = 'Load More Results';
        loadMoreBtn.onclick = () => this.loadMoreAndDisplay(this.requestToken);
        resultsContainer.appendChild(loadMoreBtn);
    }


    async displayResults(requestToken = this.requestToken) {
        const resultsContainer = document.getElementById('results-container');
        resultsContainer.innerHTML = this.loadingHtml;

        const loadState = this.loaded === null ? { ok: false, stale: false } : await this.loaded;
        if (requestToken !== this.requestToken || loadState.stale) {
            return;
        }

        if (!loadState.ok) {
            resultsContainer.innerHTML = '<div class="no-results">Error occured during compounds loading</div>';
            return;
        }

        if (this.results.length === 0) {
            resultsContainer.innerHTML = '<div class="no-results">No compounds found</div>';
            return;
        }

        resultsContainer.innerHTML = '';
        this.results.forEach((cid) => {
            resultsContainer.appendChild(this.#createResultItem(cid));
        });

        this.syncSelectedLevelInputs();

        if (!this.isEnd) {
            this.#addLoadMoreButton();
        }
    }


    showPopupNode(cid) {
        const compound = this.compounds.get(cid);
        const name = compound.name;

        let html = '<div class="compound-structure">';
        html += `<img src="${escapeAttribute(this.compounds.getSvgPath(cid))}" alt="${escapeAttribute(name)}">`;
        html += '</div>';

        if (compound.nfpa || compound.pictograms) {
            html += '<div class="compound-hazards">';
            if (compound.nfpa) {
                const createNfpaDiamond = (health, flammability, instability) => `
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
                </div>`;

                const nfpa = compound.nfpa;
                const health = 'health' in nfpa ? escapeHtml(nfpa.health) : "";
                const flammability = 'flammability' in nfpa ? escapeHtml(nfpa.flammability) : "";
                const instability = 'instability' in nfpa ? escapeHtml(nfpa.instability) : "";
                html += createNfpaDiamond(health, flammability, instability);
            }

            if (compound.pictograms) {
                html += '<div class="ghs-pictograms">';
                compound.pictograms.forEach((ghs) => {
                    html += `<img src="/assets/ghs_pictograms/${escapeAttribute(ghs)}.svg" alt="${escapeAttribute(ghs)} pictogram">`;
                });
                html += '</div>';
            }

            html += '</div>';
        }

        html += '<div class="compound-info compound-info-general">';
        const properties = Array.isArray(compound.properties) ? compound.properties : [];
        properties.forEach((entry) => {
            const label = String(entry.property_name ?? "");
            const rawValue = String(entry.property_value ?? "");
            const labelLower = label.toLowerCase();

            html += '<div class="compound-info-row">';
            html += `<span class="compound-info-label table-info">${escapeHtml(label)}</span>`;
            html += '<span class="compound-info-value table-info">';

            if (labelLower.includes('pubchem') && labelLower.includes('cid')) {
                const cidValue = rawValue.trim();
                html += `<a href="https://pubchem.ncbi.nlm.nih.gov/compound/${encodeURIComponent(cidValue)}" target="_blank" rel="noopener noreferrer">${escapeHtml(cidValue)}</a>`;
            } else if (labelLower.includes('wikipedia')) {
                const safeUrl = getSafeExternalUrl(rawValue);
                if (safeUrl) {
                    const linkLabel = safeUrl.split('/').at(-1).replace(/_/g, ' ');
                    html += `<a href="${escapeAttribute(safeUrl)}" target="_blank" rel="noopener noreferrer">${escapeHtml(linkLabel)}</a>`;
                } else {
                    html += formatMultilineHtml(rawValue);
                }
            } else {
                html += formatMultilineHtml(rawValue);
            }

            html += '</span>';
            html += '</div>';
        });
        html += '</div>';

        if (compound.description) {
            html += '<div class="compound-info compound-info-description">';
            html += `<span class="compound-info-value description">${formatParagraphsHtml(compound.description)}</span>`;
            html += '</div>';
        }

        this.popup.showPopup(html, name, 'popup-node');
    }
}
