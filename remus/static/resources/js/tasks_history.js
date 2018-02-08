$(document).ready(function () {
    if (Cookies.get("remus_history") === undefined) {
        Cookies.set("remus_history", []);
    }
    setItemsFromHistory();
});

function archivePerformedTask(record) {
    addHistoryRecord(record);
    setItemsFromHistory();
}

function createHistoryRecord() {
    return {
        "date": moment().format("HH:mm:ss MM/DD/YYYY"),
        "genome": $('#genome_select').val(),
        "genes": $('#genes_select').val(),
        "tissues": $('#tissues_select').val(),
        "operation": $('#operation_select').val()
    }
}

function objToString(obj) {
    var str = '';
    for (var p in obj) {
        if (obj.hasOwnProperty(p)) {
            str += '  ' + p + '\t:\t' + obj[p] + '</br>';
        }
    }
    return str;
}

function setItemsFromHistory() {
    $("#history_list").empty();
    $.each(Cookies.getJSON("remus_history").reverse(), function (index, record) {
        var date_button = "<button class='record-button' onclick='performFromHistory(this)'>" + record.date + "</button>";
        var delete_button = "<button onclick='removeFromHistory(this)'>x</button>";
        var details_popup = "<div class='popup arrow_box' style='display: none'>" + objToString(record) + "</div>";
        var list_element = "<li>" + date_button + delete_button + details_popup + "</li>";
        $("#history_list").append(list_element);
        $(document).on("mouseover", ".record-button", function (e) {
            var current_li = $(e.target);
            var x = current_li.parent().children(".popup");
            x.stop(true).fadeIn(600);
        });

        $(document).on("mouseover", ".popup", function (e) {
            var current_li = $(e.target);
            var x = current_li.parent().children(".popup");
            x.stop(true).fadeIn(600);
        });

        $(document).on("mouseout", ".record-button", function (e) {
            var current_li = $(e.target);
            var x = current_li.parent().children(".popup");
            x.stop(true).fadeOut(600);
        });

        $(document).on("mouseout", ".popup", function (e) {
            var current_li = $(e.target);
            var x = current_li.parent().children(".popup");
            x.stop(true).fadeOut(600);
        });
    });
}

function performFromHistory(element) {
    $.each(Cookies.getJSON("remus_history").reverse(), function (index, record) {
        if (element.textContent === record.date) {
            delete record.date;
            perform_operation(serialize(record))
        }
    });
}

function serialize(obj) {
    var str = [];
    for (var p in obj) {
        if (obj.hasOwnProperty(p)) {
            str.push(encodeURIComponent(p) + "=" + encodeURIComponent(obj[p]));
        }
    }
    return str.join("&");
}

function removeFromHistory(element) {
    var list_row = $(element).parent();
    var date = list_row.children(".record-button")[0].textContent;
    delHistoryRecordByDate(date);
    list_row.remove()
}

function delHistoryRecordByDate(date_to_remove) {
    var date = date_to_remove;
    var params_history = Cookies.getJSON("remus_history");
    $.each(params_history, function (index, record) {
        if (record !== undefined && date === record.date) {
            params_history.splice(index, 1)
        }
    });
    Cookies.set("remus_history", params_history);
}

function addHistoryRecord(record) {
    var params_history = Cookies.getJSON("remus_history").reverse();
    if (params_history.length === 0 || record.date !== params_history[0]["date"]) {
        params_history.push(record);
        Cookies.set("remus_history", params_history);
    }
}
