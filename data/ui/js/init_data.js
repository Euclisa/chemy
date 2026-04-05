let originalMainContent = '';
let compounds = new Compounds();
let reactions = new Reactions();
let popup = new Popup();
let catalog = new Catalog(compounds, popup);
let graph = new Graph(catalog, reactions);



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
        graph.max_paths = kValue.textContent === 'All' ? 99999 : Number(kValue.textContent);
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
