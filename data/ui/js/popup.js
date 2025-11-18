class Popup {
    constructor() {}

    showPopup(htmlContent, titleStr, subClass) {
        const popup = document.createElement('div');
        popup.classList.add('popup', subClass);

        popup.innerHTML = `
            <div class="popup-header">
                <span id="popup-title">Title</span>
                <span class="popup-close" onclick="this.closest('.popup').remove()">x</span>
            </div>
            <div class="popup-content" id="popup-content">
            Content
            </div>
        `;

        document.body.appendChild(popup);

        const title = popup.querySelector('#popup-title');
        const content = popup.querySelector('#popup-content');

        content.innerHTML = htmlContent;
        title.textContent = titleStr;

        popup.style.display = 'block';
        popup.style.position = 'absolute';
        popup.style.visibility = 'visible';
        popup.style.opacity = '1';

        const rect = popup.getBoundingClientRect();
        popup.style.left = (window.innerWidth / 2 - rect.width / 2) + 'px';
        popup.style.top = (window.innerHeight / 2 - rect.height / 2) + 'px';

        this.#makeDraggable(popup);
    }


     #makeDraggable(element) {
        const header = element.querySelector('.popup-header');
        header.addEventListener('mousedown', function(e) {
            e.preventDefault();
            let shiftX = e.clientX - element.getBoundingClientRect().left;
            let shiftY = e.clientY - element.getBoundingClientRect().top;

            function moveAt(pageX, pageY) {
                element.style.left = pageX - shiftX + 'px';
                element.style.top = pageY - shiftY + 'px';
            }

            function onMouseMove(e) {
                moveAt(e.pageX, e.pageY);
            }

            document.addEventListener('mousemove', onMouseMove);

            document.addEventListener('mouseup', function() {
                document.removeEventListener('mousemove', onMouseMove);
            }, {once: true});
        });
    }
}
