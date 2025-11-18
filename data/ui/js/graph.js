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
    }

    async submit() {
        showLoading('Computing graph...');

        this.edgeToRids.clear();

        const sources = Array.from(this.catalog.sourceCids);
        const targets = Array.from(this.catalog.targetCids);

        const post_body = {
            "sources": sources,
            "targets": targets,
            "max_paths": this.max_paths,
            "max_cost": this.max_cost,
            "primary_only": this.primary_only
        };

        const response = await fetch("/api/build_graph", {
            method: "POST",
            headers: {
                "Content-Type": "application/json"
            },
            body: JSON.stringify(post_body)
        });

        if (!response.ok) throw new Error(`HTTP error during graph fetch! status: ${response.status}`);

        const data = await response.json();

        let allCids = new Set();
        let allRids = new Set();
        let nodes = [];
        let links = [];
        const targetMap = new Map();

        console.log(data);
        const graph = data["graph"];

        graph.forEach(d => {
            const cid = d.cid;
            nodes.push({
                cid,
                primary: d.primary
            });
            targetMap.set(cid, new Set(d.targets));
        });

        graph.forEach(d => {
            const sourceCid = d.cid;
            d.targets.forEach(target => {
                const targetCid = target.cid;
                const edge = this.getEdgeId(sourceCid, targetCid);
                if (!this.edgeToRids.has(edge)) {
                    const isBi = targetMap.get(targetCid)?.has(sourceCid);
                    links.push({ source: sourceCid, target: targetCid, type: isBi ? "bi" : "directed", primary: d.primary});
                    
                    const currRids = new Set(target.reactions);
                    currRids.forEach(rid => allRids.add(rid));
                    this.edgeToRids.set(edge, currRids);

                    allCids.add(targetCid);
                }
                else {
                    const rids = this.edgeToRids.get(edge);
                    for (const r of target.reactions)
                        rids.add(r);
                }
            });
            allCids.add(sourceCid);
        });
        
        if(!await this.reactions.loadReactions(Array.from(allRids)))
            throw new Error("Failed to load reactions needed for graph constructions");

        for(const rid of allRids) {
            const reaction = this.reactions.get(rid);
            reaction["reactants"].forEach(entry => allCids.add(entry.cid));
            reaction["products"].forEach(entry => allCids.add(entry.cid));
        }


        if(!await this.catalog.compounds.loadCompounds(Array.from(allCids)))
            throw new Error("Failed to load compounds needed for graph constructions");

        nodes.forEach(node => {
            const cid = node.cid;
            const comp = this.catalog.compounds.get(cid);

            const maxNameLen = 50;
            const fullName = comp ? comp.name : 'Unknown';
            node.name = fullName.length > maxNameLen
                ? fullName.slice(0, maxNameLen) + '...'
                : fullName;

            node.organic = comp ? comp.organic : false;
        });


        this.#renderGraph(nodes, links);

        hideLoading();
    }

    getEdgeId(sourceCid, targetCid) {
        return [sourceCid, targetCid].sort().join('|');
    }

    #renderGraph(nodes, links) {
        // 1. Stop old simulation
        if (this.simulation) {
            this.simulation.stop();
        }

        const main = d3.select("#main");
        main.html("");  // safe now

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

        // 3. Create tooltip only ONCE
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
            .force("link", d3.forceLink(links).id(d => d.cid).distance(100))
            .force("charge", d3.forceManyBody().strength(-200))
            .force("center", d3.forceCenter(width / 2, height / 2));

        const link = g.append("g")
            .attr("stroke", "#b8a8d9")
            .selectAll("line")
            .data(links)
            .join("line")
            .attr("stroke-width", 2)
            .attr("marker-end", d => "url(#arrow)")
            .attr("marker-start", d => d.type === 'bi' ? "url(#arrow-start)" : null)
            .attr("stroke-opacity", d => d.primary ? 0.7 : 0.3)
            .style("cursor", "pointer")
            .on("mouseover", function(event, d) {
                d3.select(this).attr("stroke-width", 4).attr("stroke", "#9d89c7");
            })
            .on("mouseout", function(event, d) {
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
            .attr("opacity", d => d.primary ? defaultOpacity : secondaryOpacity )
            .style("cursor", "pointer")
            .call(d3.drag()
                .on("start", (event, d) => {
                    if (!event.active) this.simulation.alphaTarget(0.3).restart();
                    d.fx = d.x;
                    d.fy = d.y;
                })
                .on("drag", (event, d) => {
                    d.fx = event.x;
                    d.fy = event.y;
                })
                .on("end", (event, d) => {
                    if (!event.active) this.simulation.alphaTarget(0);
                    d.fx = null;
                    d.fy = null;
                })
            )
            .on("mouseover", (event, d) => {
                const element = d3.select(event.currentTarget);

                element
                    .attr('opacity', defaultOpacity)
                    .select("circle, rect, polygon")
                    .attr("r", 8).attr("width", 16)
                    .attr("height", 16)
                    .attr("points", "0,-8 8,8 -8,8");
                    hoverTimeout = setTimeout(() => {
                        this.tooltip.transition()
                            .duration(200)
                            .style("opacity", .9);
                        this.tooltip.html(`<img src="${this.catalog.compounds.getSvgPath(d.cid)}" alt="${d.name}"><div class="tooltip-name">${d.name}</div>`)
                            .style("left", (event.pageX + 10) + "px")
                            .style("top", (event.pageY - 28) + "px");
                    }, 500);
            })
            .on("mousemove", (event, d) => {
                this.tooltip.style("left", (event.pageX + 10) + "px")
                    .style("top", (event.pageY - 28) + "px");
            })
            .on("mouseout", (event, d) => {
                const element = d3.select(event.currentTarget);
                
                element
                    .attr('opacity', d => d.primary ? defaultOpacity : secondaryOpacity)
                    .select("circle, rect, polygon")
                    .attr("r", 5)
                    .attr("width", 10)
                    .attr("height", 10)
                    .attr("points", "0,-5 5,5 -5,5");
                    clearTimeout(hoverTimeout);
                    this.tooltip.transition()
                        .duration(500)
                        .style("opacity", 0);
            })
            .on("click", (event, d) => {
                this.catalog.showPopupNode(d.cid);
            });
        
        const catalog = this.catalog;

        nodeGroups.each(function(d) {
            const sel = d3.select(this);
            const color = d.organic ? "#7bc67b" : "#6ba3d6";
            if (catalog.sourceCids.has(d.cid)) {
                sel.append("rect")
                    .attr("x", -5)
                    .attr("y", -5)
                    .attr("width", 10)
                    .attr("height", 10)
                    .attr("fill", color);
            } else if (catalog.targetCids.has(d.cid)) {
                sel.append("polygon")
                    .attr("points", "0,-5 5,5 -5,5")
                    .attr("fill", color);
            } else {
                sel.append("circle")
                    .attr("r", 5)
                    .attr("fill", color);
            }
        });

        const labels = nodeGroups.append("text")
            .text(d => d.name)
            .attr("font-size", 10)
            .attr("text-anchor", "middle")
            .attr("dy", 20)
            .attr("fill", "#4a4a6a")
            .attr("font-weight", 500);

        this.simulation.on("tick", () => {
            link
                .attr("x1", d => {
                    const dx = d.target.x - d.source.x;
                    const dy = d.target.y - d.source.y;
                    const dr = Math.sqrt(dx * dx + dy * dy);
                    const normX = dx / dr;
                    const sourcePadding = d.type === 'bi' ? 8 : 0;
                    return d.source.x + (normX * sourcePadding);
                })
                .attr("y1", d => {
                    const dx = d.target.x - d.source.x;
                    const dy = d.target.y - d.source.y;
                    const dr = Math.sqrt(dx * dx + dy * dy);
                    const normY = dy / dr;
                    const sourcePadding = d.type === 'bi' ? 8 : 0;
                    return d.source.y + (normY * sourcePadding);
                })
                .attr("x2", d => {
                    const dx = d.target.x - d.source.x;
                    const dy = d.target.y - d.source.y;
                    const dr = Math.sqrt(dx * dx + dy * dy);
                    const normX = dx / dr;
                    const targetPadding = (d.type === 'directed' || d.type === 'bi') ? 8 : 0;
                    return d.target.x - (normX * targetPadding);
                })
                .attr("y2", d => {
                    const dx = d.target.x - d.source.x;
                    const dy = d.target.y - d.source.y;
                    const dr = Math.sqrt(dx * dx + dy * dy);
                    const normY = dy / dr;
                    const targetPadding = (d.type === 'directed' || d.type === 'bi') ? 8 : 0;
                    return d.target.y - (normY * targetPadding);
                });

            nodeGroups
                .attr("transform", d => `translate(${d.x}, ${d.y})`);
        });

        const zoom = d3.zoom()
            .scaleExtent([0.5, 4])
            .on("zoom", ({transform}) => {
                g.attr("transform", transform);
            });

        svg.call(zoom);
    }


    showPopupEdge(edge_data) {
        const sourceCid = edge_data.source.cid;
        const targetCid = edge_data.target.cid;
        const sourceComp = compounds.get(sourceCid);
        const targetComp = compounds.get(targetCid);
        const sourceName = sourceComp.name;
        const targetName = targetComp.name;
        
        const edgeId = this.getEdgeId(sourceCid, targetCid);
        const reverseEdgeId = this.getEdgeId(targetCid, sourceCid);
        const forwardReactionIDs = this.edgeToRids.get(edgeId) || [];
        const reverseReactionIDs = edge_data.type === 'bi' ? (this.edgeToRids.get(reverseEdgeId) || []) : [];
        
        let html = '';

        const generateReactionParticipants = (participants, balanced) => {
            let reactants_htmls = [];
            for (const entry of participants) {
                const cid = entry.cid;
                const compound = this.catalog.compounds.get(cid);
                const name = compound.name;
                const svgPath = this.catalog.compounds.getSvgPath(cid);
                let curr_html = ""
                curr_html += `<div class="reaction-participant">`;
                curr_html += `<img class="reaction-svg-image" src="${svgPath}" alt="${name}" data-cid="${cid}">`;
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

            let item_html = "";
            const reaction = this.reactions.get(rid);
            if (reaction) {
                const balanced = reaction.balanced;
                item_html += '<div class="reaction-item">';
                item_html += `<div class="reaction-equation">`;
                item_html += generateReactionParticipants(reaction.reactants, balanced);
                item_html += '<span class="reaction-sep reaction-sep-arrow">→</span>';
                item_html += generateReactionParticipants(reaction.products, balanced);
                item_html += '</div>';

                const description = reaction.description;
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
                    item_html += `<img src="assets/ord.svg" class="reaction-ord-label" title="Sourced from Open Reaction Database">`
                }

                item_html += '</div>';
            }

            return item_html;
        };

        const generateReactionList = (rids) => {
            let reactionsHtml = "";
            if (rids.size > 0) {
                rids.forEach(rid => {
                    reactionsHtml += generateReactionItem(rid);
                });
            } else {
                reactionsHtml += '<div class="reaction-item">No reactions found</div>';
            }

            return reactionsHtml;
        };
        
        if (edge_data.type === 'bi') {
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

        this.catalog.popup.showPopup(html, 'Reactions', 'popup-edge');

        document.querySelectorAll('.reaction-svg-image').forEach(el => {
            el.onclick = () => {
                const cid = Number(el.dataset.cid);
                this.#showCompoundInfoPopup(cid);
            };
        });
    }


    #showCompoundInfoPopup(cid) {
        const name = this.catalog.compounds.get(cid).name;
        this.catalog.showPopupNode(cid);
    }
}
