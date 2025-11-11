class Compounds {
    constructor() {
        this.compoundsMap = new Map();
        this.apiEndpoint = "/api/compound_info";
    }

    async loadCompounds(cids) {
        const cidsToFetch = cids.filter(cid => !this.compoundsMap.has(cid));
        if (cidsToFetch.length === 0) return true;

        try {
            const response = await fetch(this.apiEndpoint, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(cidsToFetch)
            });

            if (!response.ok) throw new Error(`HTTP error! status: ${response.status}`);

            const data = await response.json();
            for (const compound of data) {
                this.compoundsMap.set(compound.cid, compound);
            }

            return true;
        } catch (error) {
            console.error("Error loading compounds:", error);
            return false;
        }
    }

    async #loadCompoundSingle(cid) {
        try {
            const response = await fetch(this.apiEndpoint, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify([cid])
            });

            if (!response.ok) throw new Error(`HTTP error! status: ${response.status}`);

            const data = await response.json();
            if(data.length != 1) throw new Error(`Expected to fetch 1 compound got ${data.length}`);

            const compound = data[0];
            this.compoundsMap.set(compound.cid, compound);
        } catch (error) {
            console.error("Error loading compounds:", error);
        }
    }

    get(cid) {
        if(!this.compoundsMap.has(cid))
            throw new Error(`Compound with CID ${cid} isn't fetched`);

        return this.compoundsMap.get(cid);
    }
}