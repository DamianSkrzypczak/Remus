

$(document).ready(function () {        
        
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
            
                if (line.trim()=="" || interval.length<3) {
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
                                            filteredRecords[j+1] = x;
                                            ++callbackCnt;
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
        
        function loadAndFilterPlainTextVcf(vcf, bedRecords) {
            
            if (vcf.name.endsWith(".gz")) {
                alert("GZipped VCF selected - please select also its Tabix index file (.tbi).");
            } else if (!vcf.name.endsWith(".vcf")) {
                alert("Selected file does not have a .vcf extension. If it is a VCF file, please rename it.");
            } else {
            
                var reader = new FileReader();
                reader.onload = function(f) {         
                    // get file content     
                    var vcf_content = f.target.result;
                    // do the filtering 
                    var filtered_vcf_content = filterPlainTextVcf(vcf_content, bedRecords);
                    // and download result
                    download(vcf.name+'.remus_filtered.vcf', filtered_vcf_content);
                };
                reader.readAsText(vcf);
            }
        }
        
        function handleFilterVcf(evt) {

            // skip if no file was selected
            var files = evt.target.files;
            if (files.length < 1 || files.length > 2) {
                alert("Please select single uncompressed VCF file, or a pair of files: BGZipped VCF and its Tabix index.");
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
                
                    loadAndFilterPlainTextVcf(files[0], bedRecords);
                
                } else if (files.length > 1) {  // filter BGZipped and Tabix indexed VCF
                    
                    loadAndFilterTabixedVcf(files, bedRecords);
                    
                }
            }
            xhr.send();
            
        }
        
        document.getElementById('filter-vcf').addEventListener('change', handleFilterVcf, false); 


});
