$(document).ready(function () {
    $('#genome_select').select2({
        width: '100%',
        minimumResultsForSearch: -1,
        placeholder: "Select genome...",
        ajax: {
            url: $SCRIPT_ROOT + "/genomes",
            data: function () {
                return {};
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
});

$(document).ready(function () {
    $('#genes_select').select2({
        width: '100%',
        minimumResultsForSearch: -1,
        placeholder: "Select genes...",
        ajax: {
            url: $SCRIPT_ROOT + "/genes",
            data: function (params) {
                if (params["term"] === undefined) {
                    return {};
                } else {
                    return {
                        "genome": $("#genome_select").find(":selected").text(),
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
});

$(document).ready(function () {
    $('#tissues_select').select2({
        width: '100%',
        minimumResultsForSearch: -1,
        placeholder: "Select tissues...",
        ajax: {
            url: $SCRIPT_ROOT + "/tissues",
            data: function () {
                return {};
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
});

$(document).ready(function () {
    $('#operation_select').select2({
        width: '100%',
        minimumResultsForSearch: -1,
        placeholder: "Select operation...",
        ajax: {
            url: $SCRIPT_ROOT + "/operations",
            data: function () {
                return {};
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
});
