class Graph {
    constructor(catalog, reactions) {
        this.catalog = catalog;
        this.reactions = reactions;
        this.max_paths = null;
        this.max_cost = null;
        this.primary_only = false;
        this.edgeToRids = new Map();

        this.simulation = null;
        this.tooltip = null;
        this.hoverTimeout = null;

        this.submitButton = document.getElementById('submit-paths-button');
        this.submitButton.addEventListener('click', this.submit.bind(this));
    }


    async submit() {
        showLoading('Computing graph...');

        try {
            this.edgeToRids.clear();

            const sources = Array.from(this.catalog.sourceCids);
            const targets = Array.from(this.catalog.targetCids);
            if (sources.length === 0 || targets.length === 0) {
                throw new Error('Select at least one source and one target compound before building a graph.');
            }

            const response = await fetch("/api/build_graph", {
                method: "POST",
                headers: {
                    "Content-Type": "application/json"
                },
                body: JSON.stringify({
                    sources,
                    targets,
                    max_paths: this.max_paths,
                    max_cost: this.max_cost,
                    primary_only: this.primary_only
                })
            });

            if (!response.ok) {
                throw new Error(`HTTP error during graph fetch! status: ${response.status}`);
            }

            const data = await response.json();
            const graphData = Array.isArray(data.graph) ? data.graph : [];
            const targetMap = new Map();
            const pairToLink = new Map();
            const nodes = [];
            const allCids = new Set();
            const allRids = new Set();

            graphData.forEach((entry) => {
                const targetsForNode = Array.isArray(entry.targets) ? entry.targets : [];
                nodes.push({
                    cid: entry.cid,
                    primary: Boolean(entry.primary)
                });
                targetMap.set(entry.cid, new Set(targetsForNode.map((target) => target.cid)));
            });

            graphData.forEach((entry) => {
                const sourceCid = entry.cid;
                const targetsForNode = Array.isArray(entry.targets) ? entry.targets : [];

                targetsForNode.forEach((target) => {
                    const targetCid = target.cid;
                    const directedEdgeId = this.getDirectedEdgeId(sourceCid, targetCid);
                    const pairEdgeId = this.getPairEdgeId(sourceCid, targetCid);
                    const reactionIds = this.edgeToRids.get(directedEdgeId) || new Set();

                    target.reactions.forEach((rid) => {
                        reactionIds.add(rid);
                        allRids.add(rid);
                    });
                    this.edgeToRids.set(directedEdgeId, reactionIds);

                    const isBi = targetMap.get(targetCid)?.has(sourceCid) ?? false;
                    const existingLink = pairToLink.get(pairEdgeId);
                    if (existingLink) {
                        existingLink.primary = existingLink.primary || Boolean(entry.primary);
                        existingLink.type = isBi ? 'bi' : existingLink.type;
                    } else {
                        const link = {
                            source: sourceCid,
                            target: targetCid,
                            type: isBi ? 'bi' : 'directed',
                            primary: Boolean(entry.primary)
                        };
                        pairToLink.set(pairEdgeId, link);
                    }

                    allCids.add(targetCid);
                });

                allCids.add(sourceCid);
            });

            if (!await this.reactions.loadReactions(Array.from(allRids))) {
                throw new Error("Failed to load reactions needed for graph construction");
            }

            for (const rid of allRids) {
                const reaction = this.reactions.get(rid);
                reaction.reactants.forEach((entry) => allCids.add(entry.cid));
                reaction.products.forEach((entry) => allCids.add(entry.cid));
            }

            if (!await this.catalog.compounds.loadCompounds(Array.from(allCids))) {
                throw new Error("Failed to load compounds needed for graph construction");
            }

            nodes.forEach((node) => {
                const compound = this.catalog.compounds.get(node.cid);
                const fullName = compound ? compound.name : 'Unknown';
                node.name = fullName.length > 50 ? `${fullName.slice(0, 50)}...` : fullName;
                node.organic = compound ? compound.organic : false;
            });

            this.#renderGraph(nodes, Array.from(pairToLink.values()));
        } catch (error) {
            console.error('Error building graph:', error);
            this.reset();
            this.catalog.popup.showPopup(
                `<p>${escapeHtml(error.message || 'Failed to build the graph.')}</p>`,
                'Graph Error',
                'popup-message'
            );
        } finally {
            hideLoading();
        }
    }


    getPairEdgeId(sourceCid, targetCid) {
        return [sourceCid, targetCid].sort().join('|');
    }


    getDirectedEdgeId(sourceCid, targetCid) {
        return `${sourceCid}|${targetCid}`;
    }


    reset() {
        if (this.simulation) {
            this.simulation.stop();
        }
        if (this.hoverTimeout) {
            clearTimeout(this.hoverTimeout);
            this.hoverTimeout = null;
        }

        this.edgeToRids.clear();
        document.getElementById('main').innerHTML = originalMainContent;
    }


    #renderGraph(nodes, links) {
        if (this.simulation) {
            this.simulation.stop();
        }

        const main = d3.select("#main");
        main.html("");

        const width = main.node().clientWidth;
        const height = main.node().clientHeight;

        const svg = main.append("svg")
            .attr("width", width)
            .attr("height", height)
            .attr("viewBox", [0, 0, width, height])
            .style("max-width", "100%")
            .style("height", "auto");

        const defs = svg.append("defs");
        defs.append("marker")
            .attr("id", "arrow")
            .attr("viewBox", "0 -3 10 6")
            .attr("refX", 10)
            .attr("refY", 0)
            .attr("markerWidth", 6)
            .attr("markerHeight", 6)
            .attr("orient", "auto")
            .append("path")
            .attr("d", "M0,-3L10,0L0,3")
            .attr("fill", "#b8a8d9");

        defs.append("marker")
            .attr("id", "arrow-start")
            .attr("viewBox", "0 -3 10 6")
            .attr("refX", 10)
            .attr("refY", 0)
            .attr("markerWidth", 6)
            .attr("markerHeight", 6)
            .attr("orient", "auto-start-reverse")
            .append("path")
            .attr("d", "M0,-3L10,0L0,3")
            .attr("fill", "#b8a8d9");

        const g = svg.append("g");
        if (!this.tooltip) {
            this.tooltip = d3.select("body").append("div")
                .attr("class", "tooltip")
                .style("opacity", 0)
                .style("position", "absolute")
                .style("pointer-events", "none")
                .style("background", "white")
                .style("border", "1px solid #ccc")
                .style("border-radius", "4px")
                .style("padding", "8px")
                .style("box-shadow", "0 2px 10px rgba(0,0,0,0.1)")
                .style("z-index", 10000);
        }

        this.simulation = d3.forceSimulation(nodes)
            .force("link", d3.forceLink(links).id((d) => d.cid).distance(100))
            .force("charge", d3.forceManyBody().strength(-200))
            .force("center", d3.forceCenter(width / 2, height / 2));

        const link = g.append("g")
            .attr("stroke", "#b8a8d9")
            .selectAll("line")
            .data(links)
            .join("line")
            .attr("stroke-width", 2)
            .attr("marker-end", "url(#arrow)")
            .attr("marker-start", (d) => d.type === 'bi' ? "url(#arrow-start)" : null)
            .attr("stroke-opacity", (d) => d.primary ? 0.7 : 0.3)
            .style("cursor", "pointer")
            .on("mouseover", function() {
                d3.select(this).attr("stroke-width", 4).attr("stroke", "#9d89c7");
            })
            .on("mouseout", function() {
                d3.select(this).attr("stroke-width", 2).attr("stroke", "#b8a8d9");
            })
            .on("click", (event, d) => {
                this.showPopupEdge(d);
            });

        const secondaryOpacity = 0.3;
        const defaultOpacity = 1.0;
        const nodeGroups = g.append("g")
            .selectAll("g")
            .data(nodes)
            .join("g")
            .attr("opacity", (d) => d.primary ? defaultOpacity : secondaryOpacity)
            .style("cursor", "pointer")
            .call(d3.drag()
                .on("start", (event, d) => {
                    if (!event.active) {
                        this.simulation.alphaTarget(0.3).restart();
                    }
                    d.fx = d.x;
                    d.fy = d.y;
                })
                .on("drag", (event, d) => {
                    d.fx = event.x;
                    d.fy = event.y;
                })
                .on("end", (event, d) => {
                    if (!event.active) {
                        this.simulation.alphaTarget(0);
                    }
                    d.fx = null;
                    d.fy = null;
                }))
            .on("mouseover", (event, d) => {
                const element = d3.select(event.currentTarget);
                element
                    .attr('opacity', defaultOpacity)
                    .select("circle, rect, polygon")
                    .attr("r", 8)
                    .attr("width", 16)
                    .attr("height", 16)
                    .attr("points", "0,-8 8,8 -8,8");

                this.hoverTimeout = setTimeout(() => {
                    this.tooltip.transition()
                        .duration(200)
                        .style("opacity", 0.9);

                    this.tooltip.html(
                        `<img src="${escapeAttribute(this.catalog.compounds.getSvgPath(d.cid))}" alt="${escapeAttribute(d.name)}"><div class="tooltip-name">${escapeHtml(d.name)}</div>`
                    )
                        .style("left", `${event.pageX + 10}px`)
                        .style("top", `${event.pageY - 28}px`);
                }, 500);
            })
            .on("mousemove", (event) => {
                this.tooltip
                    .style("left", `${event.pageX + 10}px`)
                    .style("top", `${event.pageY - 28}px`);
            })
            .on("mouseout", (event, d) => {
                const element = d3.select(event.currentTarget);
                element
                    .attr('opacity', d.primary ? defaultOpacity : secondaryOpacity)
                    .select("circle, rect, polygon")
                    .attr("r", 5)
                    .attr("width", 10)
                    .attr("height", 10)
                    .attr("points", "0,-5 5,5 -5,5");

                clearTimeout(this.hoverTimeout);
                this.hoverTimeout = null;
                this.tooltip.transition()
                    .duration(500)
                    .style("opacity", 0);
            })
            .on("click", (event, d) => {
                this.catalog.showPopupNode(d.cid);
            });

        const catalog = this.catalog;
        nodeGroups.each(function(d) {
            const selection = d3.select(this);
            const color = d.organic ? "#7bc67b" : "#6ba3d6";

            if (catalog.sourceCids.has(d.cid)) {
                selection.append("rect")
                    .attr("x", -5)
                    .attr("y", -5)
                    .attr("width", 10)
                    .attr("height", 10)
                    .attr("fill", color);
            } else if (catalog.targetCids.has(d.cid)) {
                selection.append("polygon")
                    .attr("points", "0,-5 5,5 -5,5")
                    .attr("fill", color);
            } else {
                selection.append("circle")
                    .attr("r", 5)
                    .attr("fill", color);
            }
        });

        nodeGroups.append("text")
            .text((d) => d.name)
            .attr("font-size", 10)
            .attr("text-anchor", "middle")
            .attr("dy", 20)
            .attr("fill", "#4a4a6a")
            .attr("font-weight", 500);

        this.simulation.on("tick", () => {
            link
                .attr("x1", (d) => {
                    const dx = d.target.x - d.source.x;
                    const dy = d.target.y - d.source.y;
                    const distance = Math.sqrt(dx * dx + dy * dy) || 1;
                    const sourcePadding = d.type === 'bi' ? 8 : 0;
                    return d.source.x + (dx / distance) * sourcePadding;
                })
                .attr("y1", (d) => {
                    const dx = d.target.x - d.source.x;
                    const dy = d.target.y - d.source.y;
                    const distance = Math.sqrt(dx * dx + dy * dy) || 1;
                    const sourcePadding = d.type === 'bi' ? 8 : 0;
                    return d.source.y + (dy / distance) * sourcePadding;
                })
                .attr("x2", (d) => {
                    const dx = d.target.x - d.source.x;
                    const dy = d.target.y - d.source.y;
                    const distance = Math.sqrt(dx * dx + dy * dy) || 1;
                    const targetPadding = (d.type === 'directed' || d.type === 'bi') ? 8 : 0;
                    return d.target.x - (dx / distance) * targetPadding;
                })
                .attr("y2", (d) => {
                    const dx = d.target.x - d.source.x;
                    const dy = d.target.y - d.source.y;
                    const distance = Math.sqrt(dx * dx + dy * dy) || 1;
                    const targetPadding = (d.type === 'directed' || d.type === 'bi') ? 8 : 0;
                    return d.target.y - (dy / distance) * targetPadding;
                });

            nodeGroups.attr("transform", (d) => `translate(${d.x}, ${d.y})`);
        });

        svg.call(d3.zoom()
            .scaleExtent([0.5, 4])
            .on("zoom", ({ transform }) => {
                g.attr("transform", transform);
            }));
    }


    showPopupEdge(edgeData) {
        const sourceCid = edgeData.source.cid;
        const targetCid = edgeData.target.cid;
        const sourceName = this.catalog.compounds.get(sourceCid).name;
        const targetName = this.catalog.compounds.get(targetCid).name;
        const forwardReactionIDs = this.edgeToRids.get(this.getDirectedEdgeId(sourceCid, targetCid)) || new Set();
        const reverseReactionIDs = edgeData.type === 'bi'
            ? (this.edgeToRids.get(this.getDirectedEdgeId(targetCid, sourceCid)) || new Set())
            : new Set();

        const generateReactionParticipants = (participants, balanced) => {
            const participantHtmls = [];

            participants.forEach((entry) => {
                const compound = this.catalog.compounds.get(entry.cid);
                const name = compound.name;
                let participantHtml = '<div class="reaction-participant">';
                participantHtml += `<img class="reaction-svg-image" src="${escapeAttribute(this.catalog.compounds.getSvgPath(entry.cid))}" alt="${escapeAttribute(name)}" data-cid="${entry.cid}">`;
                participantHtml += `<span class="reaction-participant-name" title="${escapeAttribute(name)}">${escapeHtml(name)}</span>`;
                participantHtml += '</div>';

                if (balanced) {
                    participantHtml = `<div class="reaction-participant-balanced"><span class="reaction-participant-coeff">${escapeHtml(entry.coeff)}</span>${participantHtml}</div>`;
                }

                participantHtmls.push(participantHtml);
            });

            return participantHtmls.join('<span class="reaction-sep reaction-sep-plus">+</span>');
        };

        const generateReactionItem = (rid) => {
            const reaction = this.reactions.get(rid);
            if (!reaction) {
                return "";
            }

            const getConfidenceClass = (confidence) => {
                if (confidence < 0.5) return 'reaction-confidence-low';
                if (confidence < 0.7) return 'reaction-confidence-medium';
                return 'reaction-confidence-high';
            };

            let itemHtml = '<div class="reaction-item">';
            itemHtml += '<div class="reaction-equation">';
            itemHtml += generateReactionParticipants(reaction.reactants, reaction.balanced);
            itemHtml += '<span class="reaction-sep reaction-sep-arrow">→</span>';
            itemHtml += generateReactionParticipants(reaction.products, reaction.balanced);
            itemHtml += '</div>';

            if (reaction.description) {
                itemHtml += '<div class="reaction-description">';
                itemHtml += '<strong>Description:</strong> ';
                itemHtml += formatMultilineHtml(reaction.description);
                itemHtml += '</div>';
            }

            if (typeof reaction.confidence === "number") {
                itemHtml += `<span class="reaction-confidence ${getConfidenceClass(reaction.confidence)}" title="This reaction is generated automatically and this is its validation score.">${escapeHtml(reaction.confidence.toFixed(2))}</span>`;
            } else if (reaction.source === 'ord') {
                itemHtml += '<img src="assets/ord.svg" class="reaction-ord-label" title="Sourced from Open Reaction Database">';
            }

            itemHtml += '</div>';
            return itemHtml;
        };

        const generateReactionList = (reactionIds) => {
            if (reactionIds.size === 0) {
                return '<div class="reaction-item">No reactions found</div>';
            }

            let reactionsHtml = "";
            reactionIds.forEach((rid) => {
                reactionsHtml += generateReactionItem(rid);
            });
            return reactionsHtml;
        };

        let html = '';
        if (edgeData.type === 'bi') {
            html += '<div class="popup-split">';
            html += `<div class="popup-direction"><h4>${escapeHtml(sourceName)} → ${escapeHtml(targetName)}</h4>${generateReactionList(forwardReactionIDs)}</div>`;
            html += `<div class="popup-direction"><h4>${escapeHtml(targetName)} → ${escapeHtml(sourceName)}</h4>${generateReactionList(reverseReactionIDs)}</div>`;
            html += '</div>';
        } else {
            html += `<div class="popup-direction"><h4>${escapeHtml(sourceName)} → ${escapeHtml(targetName)}</h4>${generateReactionList(forwardReactionIDs)}</div>`;
        }

        const popup = this.catalog.popup.showPopup(html, 'Reactions', 'popup-edge');
        popup.querySelectorAll('.reaction-svg-image').forEach((image) => {
            image.addEventListener('click', () => {
                this.catalog.showPopupNode(Number(image.dataset.cid));
            });
        });
    }
}
