let chemsData = [];
let edgesData = [];
let reactionsData = [];
let fuse = null;
const RESULTS_PER_PAGE = 200;
let selectedCIDs = new Set();
let cidToCompound = new Map();
let originalMainContent = '';
let hoverTimeout;
let selectedLevel = 1;
let graph_reverse = new Map()
let edgeToReactionID = new Map();
let RIDToReaction = new Map();
let currentEdges = new Set();
let secondaryNodes = new Set();
let directedEdges = new Set();
let backgroundCids;
let cidToDescription = new Map();
let ridToDescription = new Map();
let commonnesSortedCids = [];
let complexitySortedCids = [];
let curiositySortedCids = [];
let currentSortingOrder;
let cidToEdges = new Map();
let cidToCompoundProperties = new Map();
let cidToHazards = new Map();

let sortingOrder = "none";
let currentPage = 0;
let compounds = new Compounds();
let catalog = new Catalog(compounds);
let graph = new Graph(catalog);

async function loadData(fileName) {
    try {
        const response = await fetch(fileName, { cache: "no-store" });
        const text = await response.text();

        let data = [];
        if (fileName.endsWith('.jsonl')) {
            data = text.split('\n').filter(line => line.trim() !== '');
            data = data.map(line => JSON.parse(line));
        }
        else if (fileName.endsWith('.json')) {
            data = JSON.parse(text);
        }

        return data;
    } catch(error) {
        console.error(`Error during fetching: ${fileName}; ${error}`);
        return [];
    }
}

async function initializeData() {
    showLoading('Loading data...');
    const resultsContainer = document.getElementById('results-container');
    originalMainContent = document.getElementById('main').innerHTML;
    
    try {
        catalog.setQuery();
        
    } catch (error) {
        console.error('Error initializing data:', error);
        resultsContainer.innerHTML = '<div class="no-results">Error loading compound data</div>';
    }
    hideLoading();

    catalog.renderPathList();

    document.getElementById('level1').addEventListener('click', () => catalog.selectLevel(1));
    document.getElementById('level2').addEventListener('click', () => catalog.selectLevel(2));

    const kSlider = document.getElementById('k-slider');
    const kValue = document.getElementById('k-value');

    const updateMaxPaths = () => {
        kValue.textContent = kSlider.value === '11' ? 'All' : kSlider.value;
        graph.max_paths = kValue.textContent === 'All' ? 99999 : Number(kValue.textContent)
    };

    updateMaxPaths();
    kSlider.addEventListener('input', updateMaxPaths);

    const nSlider = document.getElementById('n-slider');
    const nValue = document.getElementById('n-value');

    const updateMaxCost = () => {
        nValue.textContent = nSlider.value;
        graph.max_cost = Number(nSlider.value);
    };
    
    updateMaxCost();
    nSlider.addEventListener('input', updateMaxCost);
}