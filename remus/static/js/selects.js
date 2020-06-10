DELAY = 100;
$(document).ready(function () {

    $("body").on('click', '.select2-results__group', function() {
        $(this).siblings().toggle();
      })

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

    function getClss(){
        var clss = {};
        $.ajax({
            dataType: "json",
            url: "/api/classes",
            async: false,
            success: function(data) {
                $.each(data, function (i, item) {
                    clss[item["name"]] = item["class"]
                })
            }
          });
        return clss
    }

    // $(".select2-results__options--nested").ready(function() {
    //     console.log("LOL");
    //     $(this).siblings().toggle();
    //   })



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


                var cells = []
                var tissues = []
                var organs = []
                var clss = getClss()
                $.each(data, function (i, item) {
                    var cls = clss[item]
                    if (cls === "cell"){
                        cells.push({"id": item, "text": item})
                    } else if (cls === "tissue") {
                        tissues.push({"id": item, "text": item})
                    } else if (cls === "organ") {
                        organs.push({"id": item, "text": item})
                    }
                });
                return {results: results = [
                    {
                        id: 'cells',
                        text: 'cells (' + cells.length + ` results) (click to show/hide)`,
                        disabled: true,
                        children: cells
                    },
                    {
                        id: 'tissues',
                        text: 'tissues (' + tissues.length + ` results) (click to show/hide)`,
                        disabled: true,
                        children: tissues
                    },
                    {
                        id: 'organs',
                        text: 'organs (' + organs.length + ` results) (click to show/hide)`,
                        disabled: true,
                        children: organs
                    }
                ]};
            }
        }
    });

    $("form").submit(function (e) {
        e.preventDefault();
        $('#results-table').hide();
        $('#results-loading').show();
        $('#download-result').attr("disabled", true);
        $('#download-excel').attr("disabled", true);
        $('#link-genomebrowser').attr("disabled", true);
        $.ajax({
            type: "POST",
            url: $SCRIPT_ROOT + "/api/perform",
            data: $("#main-form").serialize() + "&genome=" + get_genome(),
            success: function (result) {
                $("#results-table").html("<h3>Summary table</h3>" + result);
                $('#results-loading').hide();
                $('#results-table').show();
                $('#download-result').attr("disabled", false);
                $('#download-excel').attr("disabled", false);
                $('#link-genomebrowser').attr("disabled", false);
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

    $("#download-excel").bind("click", "doubleclick", (function (e) {
        e.preventDefault();
        window.location.href = $SCRIPT_ROOT + "/api/download_last_excel";
        return false
    }));

    function get_last_result_id() {
        var id = "";
        $.ajax({
            dataType: "text",
            url: "/api/last_result_id",
            async: false,
            success: function(data) {
                id=data
            }
          });
        return id
    }

    function get_first_gene() {
       return $('#select-genes').select2("val")[0];
    }

    $("#link-genomebrowser").bind("click", "doubleclick", (function (e) {
        e.preventDefault();

        var gene = get_first_gene()
        if (gene === void(0)) {
            gene="chr1"
        }

        // http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=GENOME_BUILD&position=GENE&hgt.customText=http://remus.btm.umed.pl/api/dowload_by_id/ID
        var url = "http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=" + get_genome() +
                               "&position=" + gene +
                               "&hgt.customText=http://" + window.location.host + "/api/download_by_id/" +
                               get_last_result_id();
        var win = window.open(url, '_blank');
        if (win) { //Browser has allowed it to be opened
            win.focus();
        } else {  //Browser has blocked it
            alert('Please allow popups for this website');
        }
        return false
    }));

    $(function () {
        $('[data-toggle="popover"]').popover()
    });




    /*
     *
     *  VCF filtering
     *
     */


        function download(filename, text) {
            var element = document.createElement('a');
            element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
            element.setAttribute('download', filename);

            element.style.display = 'none';
            document.body.appendChild(element);

            element.click();

            document.body.removeChild(element);
        }


        function filterPlainTextVcf(input_vcf_content, intervals) {

            //console.log('in_filter:'+input_vcf_content);
            //console.log('in_filter:'+bed_content);

            var trimChr = function(chr) { return new String(chr).toLowerCase().replace('chr',''); }

            var filtered_vcf_content="";
            var vcfLines = input_vcf_content.split(/\r\n|\n/);

            var vcf_index=0;
            var interval_index=0;
            while ( interval_index < intervals.length && vcf_index < vcfLines.length) {

                line = vcfLines[vcf_index];
                interval = intervals[interval_index];

                // if (vcf_index % 100 == 0) {  console.log(vcf_index); }

                if (line.trim()=="") {
                    // end of file
                    break;
                }

                if (interval.length<3) {
                    alert("Incorrect interval: ", interval);
                    break;
                }
                i_chr = parseInt(trimChr(interval[0]));
                i_start = parseInt(interval[1]);
                i_end = parseInt(interval[2]);

                if (line.startsWith("#")) {
                    filtered_vcf_content += line+"\n";
                    vcf_index++;
                } else {
                    vcf_records = line.split('\t');
                    pos = parseInt(vcf_records[1]);
                    vcf_chr = parseInt(trimChr(vcf_records[0]));

                    if (vcf_chr < i_chr ||
                        (vcf_chr == i_chr && pos < i_start)) {
                        // vcf record < interval    =>    move to the next record
                        vcf_index++;

                    } else if ((vcf_chr == i_chr && pos >= i_end) ||
                               vcf_chr > i_chr) {
                        // vcf_record > interval    =>    move to the next interval
                        interval_index++;
                    } else {
                        // console.log(line + "\n" + interval + "\n" + vcf_chr+" "+pos + "\n" + i_chr+" "+i_start+" "+i_end +"\n" + (pos < i_end) + " " + (pos>=i_start))

                        // variant in interval boundary
                        filtered_vcf_content += line+"\n";
                        vcf_index++;
                    }
                }
            }

            return filtered_vcf_content;
        }



        function makeVcfString(list) {
            var result = [];
            list.forEach(function(e) {
                result = result.concat(e);
            });
            return result.join('\n');
        }

        function filterTabixedVcfCallbackMerger(reader, bed) {

            var filteredVcfName = reader.theFile.name + ".remus_filtered.vcf"
            var filteredRecords = [];
            var callbackCnt = 0;

            reader.getHeader(function(x) {
                                filteredRecords[0] = x;
                                ++callbackCnt;
                                // console.log("Got header. Cbckcnt:", callbackCnt);
                                if (callbackCnt == bed.length+1) {
                                    download(filteredVcfName, makeVcfString(filteredRecords));
                                }
                            });

            for (i=0; i<bed.length; i++) {
                reg = bed[i];
                reader.getRecords(reg[0], reg[1], reg[2],
                                  (function () {
                                        var j = i;
                                        return function(x) {
                                            if (x.length > 0) {
                                                console.log('x=[',x,']');
                                            }
                                            filteredRecords[j+1] = x;
                                            ++callbackCnt;
                                            //console.log("Got ", j, ". Cbckcnt:", callbackCnt);
                                            if (callbackCnt == bed.length+1) {
                                                download(filteredVcfName, makeVcfString(filteredRecords));
                                            }
                                        }
                                   })() );
            }
        }


        function loadAndFilterTabixedVcf(files, bedRecords) {

            var isTabi = function(f) {return f.name.endsWith(".tbi");};
            var tabixFile = isTabi(files[0]) ? files[0] : files[1];
            var vcfFile = isTabi(files[1]) ? files[0] : files[1];

            if (tabixFile == vcfFile || !vcfFile.name.endsWith(".gz")) {
                alert("Please select single uncompressed VCF file, or a pair of files: BGZipped VCF and its Tabix index.");
            } else {
                // console.log("Tabix: " +tabixFile.name + "\nVCF: " + vcfFile.name);
                vcfR = new readBinaryVCF(tabixFile, vcfFile,
                                    (function() {
                                        return function(vcfReader) {
                                                    return filterTabixedVcfCallbackMerger(vcfReader, bedRecords);
                                                }
                                    })());
            }
        }

        function loadAndFilterPlainTextVcf(vcf, bedRecords, call_when_done) {

            if (vcf.name.endsWith(".gz")) {
                //alert("GZipped VCF selected - please select also its Tabix index file (.tbi).");
                alert("Seems like you selected GZipped VCF file. Please decompress it and try again.");
            } else if (!vcf.name.endsWith(".vcf")) {
                alert("Selected file does not have a .vcf extension. If it is a VCF file, please rename it.");
            } else {

                var reader = new FileReader();
                reader.onload = function(f) {
                    try {
                        // get file content
                        var vcf_content = f.target.result;
                        // do the filtering
                        var filtered_vcf_content = filterPlainTextVcf(vcf_content, bedRecords);
                        // and download result
                        download(vcf.name+'.remus_filtered.vcf', filtered_vcf_content);
                    } catch (e) {
                        if (e.message == 'allocation size overflow') {
                            alert('The VCF file is too large. Divide it in parts or try a different browser');
                        } else {
                            alert(e);
                        }
                    } finally {
                        call_when_done();
                    }
                };
                reader.readAsText(vcf);
            }
        }

        function handleFilterVcf(evt, call_when_done) {

            // skip if no file was selected
            var files = evt.target.files;
            if (files.length < 1 || files.length > 2) {
                //alert("Please select single uncompressed VCF file, or a pair of files: BGZipped VCF and its Tabix index.");
                alert("Please select single uncompressed VCF file.");
                return;
            }

            // download the result BED first (while user selects VCF files)
            var xhr = new XMLHttpRequest();
            xhr.open("GET", $SCRIPT_ROOT + "/api/download_last");
            xhr.responseType = "text";
            xhr.onload = function() {

                var bed = xhr.response;
                //console.log('bed1: ' + bed);

                // parse bed file into matrix of coordinates
                bedRecords = [];
                bed.split(/\r\n|\n/).forEach(function(x) {
                                                r = x.split("\t")
                                                if (r.length>=3) {
                                                    bedRecords.push(r)
                                                }
                                            });
                //console.log(bedRecords);


                // load and filter user's VCF
                if (files.length == 1) {   // single not-compressed VCF

                    loadAndFilterPlainTextVcf(files[0], bedRecords, call_when_done);

                } else if (files.length > 1) {  // filter BGZipped and Tabix indexed VCF

                    alert("Filtering BGZipped VCF files is currently performing suboptimally. Please use plain-text VCF file instead.");
                    // loadAndFilterTabixedVcf(files, bedRecords);

                }
            }
            xhr.send();

        }


    $("#filter-vcf").bind("change", (function (e) {
        e.preventDefault();
        $('#results-table').hide();
        $('#results-loading2').show();
        handleFilterVcf(e, (function() {
            $('#results-loading2').hide();
            $('#results-table').show();
        }));

        return false
    }));



});


