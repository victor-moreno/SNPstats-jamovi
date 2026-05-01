'use strict';

module.exports = {

    weightsPath_creating: function(ui) {
        console.log('weightsPath fired');

        let $el = ui.weightsPath.$el;

        let $input = $el.find('input');

        console.log('input length:', $input.length);

        setTimeout(() => {
            let $input = $el.find('input');

            let $btn = $('<button type="button">Browse…</button>');
            $input.after($btn);

            $btn.on('click', () => {
                alert('clicked');
            });
        }, 0);
    }

};