
const HTML_ESCAPE_MAP = {
    '&': '&amp;',
    '<': '&lt;',
    '>': '&gt;',
    '"': '&quot;',
    "'": '&#39;'
};


function escapeHtml(value = "") {
    return String(value).replace(/[&<>"']/g, (char) => HTML_ESCAPE_MAP[char]);
}


function escapeAttribute(value = "") {
    return escapeHtml(value);
}


function formatMultilineHtml(value = "") {
    return escapeHtml(value).replace(/\n/g, '<br>');
}


function formatParagraphsHtml(value = "") {
    const paragraphs = String(value)
        .trim()
        .split(/\n\s*\n/)
        .map((paragraph) => paragraph.trim())
        .filter(Boolean);

    return paragraphs.map((paragraph) => `<p>${formatMultilineHtml(paragraph)}</p>`).join('');
}


function getSafeExternalUrl(url) {
    try {
        const parsedUrl = new URL(String(url), window.location.origin);
        if (parsedUrl.protocol === 'http:' || parsedUrl.protocol === 'https:') {
            return parsedUrl.href;
        }
    } catch (error) {
        return null;
    }

    return null;
}


function returnToBlank() {
    graph.reset();
    catalog.resetPathLists();
}


function loadStructureSVG(cid, container) {
    const svgPath = `/assets/structures/${cid}.svg`;

    fetch(svgPath)
        .then((response) => {
            if (!response.ok) {
                throw new Error(`HTTP ${response.status}`);
            }
            return response.text();
        })
        .then((svgContent) => {
            container.innerHTML = svgContent;
            container.classList.remove('loading');
        })
        .catch(() => {
            container.innerHTML = '<i class="fa fa-flask" style="color: #ccc;"></i>';
            container.classList.remove('loading');
        });
}


function showLoading(message, showCancel = false) {
    const overlay = document.getElementById('loading-overlay');
    document.getElementById('loading-message').textContent = message;
    document.getElementById('cancel-compute').style.display = showCancel ? 'block' : 'none';
    overlay.style.display = 'flex';
}

function hideLoading() {
    document.getElementById('loading-overlay').style.display = 'none';
}
