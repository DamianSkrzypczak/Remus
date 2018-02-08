$(document).ready(function () {
    $('#operation_form').on('keyup keypress', function (e) {
        var keyCode = e.keyCode || e.which;
        if (keyCode === 13) {
            e.preventDefault();
            return false;
        }
    });
});

$(document).ready(function () {
    $('#operation_form').submit(function (e) {
        perform_operation($('#operation_form').serialize());
        e.preventDefault();
        return false
    });
});

function perform_operation(json_string) {
    console.log(json_string);
    var genome = $('#genome_select').val();
    var genes = $('#genes_select').val();
    var tissues = $('#tissues_select').val();
    var operation = $('#operation_select').val();
    if (((genome !== "" && genes.length >= 1) || tissues.length >= 1) && operation !== "") {
        $("#summary").html('<img src="static/loading2.gif"/>');
        $.ajax({
            type: 'POST',
            url: $SCRIPT_ROOT + '/perform',
            data: json_string,
            success: function (data) {
                $("#summary").html(data);
                var record = createHistoryRecord();
                archivePerformedTask(record);
            }
        });
    } else {
        alert("Please choose operation and at least one gene/tissue")
    }
}



