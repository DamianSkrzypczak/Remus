DELAY = 100;
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
        // allowClear: true,
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
        minimumInputLength: -1,
        placeholder: "Select organs/tissues/cell types",
        // closeOnSelect: false,
        // allowClear: true,
        ajax: {
            delay: DELAY,
            url: $SCRIPT_ROOT + "/api/tissues",
            data: function (params) {
                if (params["term"] === undefined) {
                    return {};
                } else {
                    return {
                        "genome": get_genome(),
                        "pattern": params.term,
                        "limit": 0
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
        $('#results-table').hide();
        $('#results-loading').show();
        $('#download-result').attr("disabled", true);
        $.ajax({
            type: "POST",
            url: $SCRIPT_ROOT + "/api/perform",
            data: $("#main-form").serialize() + "&genome=" + get_genome(),
            success: function (result) {
                $("#results-table").html("<h3>Summary table</h3>" + result);
                $('#results-loading').hide();
                $('#results-table').show();
                $('#download-result').attr("disabled", false);
                $('#filter-vcf-label').removeClass('disabled');
                $('#filter-vcf').attr("disabled", false);
                $('#filter-vcf').val('');
            }
            // error: function (result) {
            //     alert('error');
            // }
        });
        return false;
    });


    $("#download-result").bind("click", "doubleclick", (function (e) {
        e.preventDefault();
        window.location.href = $SCRIPT_ROOT + "/api/download_last";
        return false
    }));


    $(function () {
        $('[data-toggle="popover"]').popover()
    });
});


