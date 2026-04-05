class Popup {
    constructor() {}

    showPopup(htmlContent, titleStr, subClass) {
        const popup = document.createElement('div');
        popup.classList.add('popup', subClass);

        popup.innerHTML = `
            <div class="popup-header">
                <span class="popup-title">Title</span>
                <button type="button" class="popup-close" aria-label="Close popup">x</button>
            </div>
            <div class="popup-content">
            Content
            </div>
        `;

        document.body.appendChild(popup);

        const title = popup.querySelector('.popup-title');
        const content = popup.querySelector('.popup-content');
        const closeButton = popup.querySelector('.popup-close');

        content.innerHTML = htmlContent;
        title.textContent = titleStr;
        closeButton.addEventListener('click', () => popup.remove());

        popup.style.display = 'block';
        popup.style.position = 'fixed';
        popup.style.visibility = 'visible';
        popup.style.opacity = '1';

        const rect = popup.getBoundingClientRect();
        popup.style.left = `${Math.max((window.innerWidth - rect.width) / 2, 16)}px`;
        popup.style.top = `${Math.max((window.innerHeight - rect.height) / 2, 16)}px`;

        this.#makeDraggable(popup);
        return popup;
    }


     #makeDraggable(element) {
        const header = element.querySelector('.popup-header');
        header.addEventListener('mousedown', function(e) {
            e.preventDefault();
            let shiftX = e.clientX - element.getBoundingClientRect().left;
            let shiftY = e.clientY - element.getBoundingClientRect().top;

            function moveAt(clientX, clientY) {
                element.style.left = `${clientX - shiftX}px`;
                element.style.top = `${clientY - shiftY}px`;
            }

            function onMouseMove(e) {
                moveAt(e.clientX, e.clientY);
            }

            document.addEventListener('mousemove', onMouseMove);

            document.addEventListener('mouseup', function() {
                document.removeEventListener('mousemove', onMouseMove);
            }, {once: true});
        });
    }
}
