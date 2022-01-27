import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/io.wdl" as IO
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS

workflow optional_unzip_file{
    # Helper workflow that returns decompressed version of file regardless of whether input file is compressed
    File input_file
    String? user_filename
    Int cpu = 1
    Int mem_gb = 1

    # Check if file is gzipped
    call UTILS.get_file_extension{
        input:
            input_file = input_file
    }

    # Unzip if nessecary
    if(get_file_extension.extension == ".gz" || get_file_extension.extension == ".gzip"){
        call IO.gunzip{
            input:
                input_file = input_file,
                user_filename = user_filename,
                cpu = cpu,
                mem_gb = mem_gb
        }
    }

    output{
        File output_file = select_first([gunzip.output_file, input_file])
    }
}