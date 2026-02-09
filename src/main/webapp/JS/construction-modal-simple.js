/**
 * ==========================================================================
 * Simple Under Construction Modal - scSAID
 * Simplified version - smaller and cleaner
 * ==========================================================================
 */

(function() {
    'use strict';

    // Configuration
    const MODAL_STORAGE_KEY = 'scSAID_construction_modal_dismissed_simple';
    const MODAL_VERSION = '1.1'; // Increment this to show modal again after updates

    /**
     * Create and inject the simplified modal HTML into the page
     */
    function createModal() {
        const modalHTML = `
            <div class="construction-modal-overlay" id="constructionModalOverlay" role="dialog" aria-labelledby="modalTitle" aria-modal="true">
                <div class="construction-modal">
                    <button class="construction-modal__close" id="modalCloseBtn" aria-label="Close modal">
                        <svg class="construction-modal__close-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <line x1="18" y1="6" x2="6" y2="18"></line>
                            <line x1="6" y1="6" x2="18" y2="18"></line>
                        </svg>
                    </button>

                    <div class="construction-modal__header">
                        <div class="construction-modal__icon-wrapper">
                            <svg class="construction-modal__icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <path d="M14.7 6.3a1 1 0 0 0 0 1.4l1.6 1.6a1 1 0 0 0 1.4 0l3.77-3.77a6 6 0 0 1-7.94 7.94l-6.91 6.91a2.12 2.12 0 0 1-3-3l6.91-6.91a6 6 0 0 1 7.94-7.94l-3.76 3.76z"></path>
                            </svg>
                        </div>
                        <h2 class="construction-modal__title" id="modalTitle">Under Construction</h2>
                    </div>

                    <div class="construction-modal__body">
                        <p class="construction-modal__message">
                            The scSAID database is currently under development. Some features may be temporarily unavailable.
                        </p>
                    </div>

                    <div class="construction-modal__footer">
                        <button class="construction-modal__button" id="modalDismissBtn">
                            Continue
                        </button>
                    </div>
                </div>
            </div>
        `;

        // Insert modal at the end of body
        document.body.insertAdjacentHTML('beforeend', modalHTML);
    }

    /**
     * Check if the modal should be displayed
     */
    function shouldShowModal() {
        try {
            const dismissedData = sessionStorage.getItem(MODAL_STORAGE_KEY);

            if (!dismissedData) {
                return true;
            }

            const parsed = JSON.parse(dismissedData);
            return parsed.version !== MODAL_VERSION;
        } catch (error) {
            // If there's any error reading storage, show the modal
            return true;
        }
    }

    /**
     * Mark the modal as dismissed
     */
    function dismissModal() {
        try {
            const dismissData = {
                version: MODAL_VERSION,
                timestamp: new Date().toISOString()
            };
            sessionStorage.setItem(MODAL_STORAGE_KEY, JSON.stringify(dismissData));
        } catch (error) {
            console.warn('Failed to save modal dismissal state:', error);
        }
    }

    /**
     * Close the modal with animation
     */
    function closeModal() {
        const overlay = document.getElementById('constructionModalOverlay');
        if (!overlay) return;

        // Add closing class for exit animation
        overlay.classList.add('is-closing');

        // Remove from DOM after animation completes
        setTimeout(() => {
            overlay.remove();
        }, 250); // Match the animation duration in CSS

        // Mark as dismissed
        dismissModal();

        // Re-enable scroll
        document.body.style.overflow = '';
    }

    /**
     * Setup event listeners
     */
    function setupEventListeners() {
        const overlay = document.getElementById('constructionModalOverlay');
        const closeBtn = document.getElementById('modalCloseBtn');
        const dismissBtn = document.getElementById('modalDismissBtn');

        if (!overlay) return;

        // Close button
        if (closeBtn) {
            closeBtn.addEventListener('click', closeModal);
        }

        // Dismiss button
        if (dismissBtn) {
            dismissBtn.addEventListener('click', closeModal);
        }

        // Close on overlay click (but not on modal content click)
        overlay.addEventListener('click', function(event) {
            if (event.target === overlay) {
                closeModal();
            }
        });

        // Close on Escape key
        document.addEventListener('keydown', function(event) {
            if (event.key === 'Escape') {
                closeModal();
            }
        });
    }

    /**
     * Initialize the modal
     */
    function init() {
        // Check if we should show the modal
        if (!shouldShowModal()) {
            return;
        }

        // Wait for DOM to be ready
        if (document.readyState === 'loading') {
            document.addEventListener('DOMContentLoaded', function() {
                createModal();
                setupEventListeners();

                // Prevent body scroll when modal is open
                document.body.style.overflow = 'hidden';
            });
        } else {
            createModal();
            setupEventListeners();

            // Prevent body scroll when modal is open
            document.body.style.overflow = 'hidden';
        }
    }

    // Initialize
    init();

})();