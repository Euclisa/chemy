class Reactions {
    constructor() {
        this.reactionsMap = new Map();
        this.apiEndpoint = "/api/reaction_info";
    }

    async loadReactions(rids) {
        const ridsToFetch = rids.filter(rids => !this.reactionsMap.has(rids));
        if (ridsToFetch.length === 0) return true;

        try {
            const response = await fetch(this.apiEndpoint, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(ridsToFetch)
            });

            if (!response.ok) throw new Error(`HTTP error! status: ${response.status}`);

            const data = await response.json();
            for (const reaction of data) {
                this.reactionsMap.set(reaction.rid, reaction);
            }

            return true;
        } catch (error) {
            console.error("Error loading reactions:", error);
            return false;
        }
    }

    get(rid) {
        if(!this.reactionsMap.has(rid))
            throw new Error(`Reaction with RID '${rid}' isn't fetched`);

        return this.reactionsMap.get(rid);
    }
}