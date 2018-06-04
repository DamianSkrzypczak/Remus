DELAY = 10;
$(document).ready(function () {

    function get_genome() {
        return window.getComputedStyle(document.querySelector('#genome-switch-label'), ':after')
            .getPropertyValue('content').replace(/['"]+/g, '')
    }

    $('#select-genes').select2({
        multiple: true,
        width: '100%',
        minimumInputLength: 1,
        placeholder: "Select genes",
        allowClear: true,
        ajax: {
            delay: DELAY,
            url: $SCRIPT_ROOT + "/api/genes",
            data: function (params) {
                if (params["term"] === undefined) {
                    return {};
                } else {
                    return {
                        "genome": get_genome(),
                        "pattern": params.term,
                        "limit": 10
                    }
                }
            },
            processResults: function (data) {
                var results = [];
                $.each(data, function (i, item) {
                    results.push({"id": item, "text": item})
                });
                return {results: results};
            }
        }
    });

    $('#select-tissues').select2({
        multiple: true,
        width: '100%',
        minimumInputLength: 0,
        placeholder: "Select tissues",
        closeOnSelect: false,
        allowClear: true,
        ajax: {
            delay: DELAY,
            url: $SCRIPT_ROOT + "/api/tissues",
            data: function (params) {
                if (params["term"] === undefined) {
                    return {};
                } else {
                    return {
                        "pattern": params.term,
                        "limit": 10
                    }
                }
            },
            processResults: function (data) {
                var results = [];
                $.each(data, function (i, item) {
                    results.push({"id": item, "text": item})
                });
                return {results: results};
            }
        }
    });

    $('#select-operation').select2({
        multiple: false,
        width: '100%',
        placeholder: "Choose operation",
        allowClear: true,
        ajax: {
            delay: DELAY,
            url: $SCRIPT_ROOT + "/api/operations",
            data: function (params) {
                if (params["term"] === undefined) {
                    return {};
                } else {
                    return {
                        "pattern": params.term,
                        "limit": 10
                    }
                }
            },
            processResults: function (data) {
                var results = [];
                $.each(data, function (i, item) {
                    results.push({"id": item, "text": item})
                });
                return {results: results};
            }
        }
    });

    $("form").submit(function (e) {
        e.preventDefault();
        var data = $("#main-form").serializeArray(); // convert form to array
        data.push({name: "genome", value: get_genome()});
        $('#results-table').hide();
        $('#results-loading').show();
        // $('#submit-all').attr("disabled", true);
        $('#download-result').attr("disabled", true);
        $.ajax({
            type: "POST",
            url: $SCRIPT_ROOT + "/api/perform",
            data: $.param(data),
            success: function (result) {
                $("#results-table").html("<p>Summary table:</p>" + result);
                $('#results-loading').hide();
                $('#results-table').show();
                // $('#submit-all').attr("disabled", false);
                $('#download-result').attr("disabled", false);
            }
            // error: function (result) {
            //     alert('error');
            // }
        });
        return false;
    });
});


