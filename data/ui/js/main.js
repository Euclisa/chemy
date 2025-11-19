(function(){
    const app = document.getElementById('app');
    const toggle = document.getElementById('toggle');
    const sidePanel = document.getElementById('sidePanel');

    let isOpen = false;

    function setOpen(open){
        if(open){
            app.classList.add('open');
            toggle.setAttribute('aria-expanded','true');
            sidePanel.setAttribute('aria-hidden','false');
        } else {
            app.classList.remove('open');
            toggle.setAttribute('aria-expanded','false');
            sidePanel.setAttribute('aria-hidden','true');
        }
    }

    setOpen(isOpen);

    toggle.addEventListener('click', function(){
        isOpen = !isOpen;
        setOpen(isOpen);
    });

    document.addEventListener('keydown', function(e){
        if(e.key === 'Escape' && isOpen){
            isOpen = false;
            setOpen(isOpen);
            toggle.focus();
        }
    });


    const clear_btn = document.getElementById('floating-btn-clear');
    clear_btn.addEventListener('click', returnToBlank);
})();

(function(){
    document.addEventListener('click', (e) => {
        const dropdown = e.target.closest('.filter-dropdown');
        
        // Close all dropdowns first
        document.querySelectorAll('.filter-dropdown').forEach(d => {
            if (d !== dropdown) d.classList.remove('active');
        });

        // Toggle current dropdown if button clicked
        if (dropdown && e.target.closest('#catalog-sort-btn')) {
            dropdown.classList.toggle('active');
        } else {
            // Clicked outside — close all
            document.querySelectorAll('.filter-dropdown').forEach(d => d.classList.remove('active'));
        }
    });
})();

document.addEventListener('DOMContentLoaded', initializeData);