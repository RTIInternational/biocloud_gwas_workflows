import "biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV

workflow test_tsv_utils{
    Array[File] tsv_inputs_to_append
    File tsv_input
    File input_filter_tsv
    File gzipped_tsv_input
    File gzipped_input_filter_tsv
    String delimiter = " "

    call TSV.tsv_append{
        input:
            tsv_inputs = tsv_inputs_to_append,
            delimiter = delimiter,
            output_filename = "append_tsv.txt",
            header = true
    }
    call TSV.tsv_filter as filter_unzipped{
        input:
            tsv_input = tsv_input,
            delimiter = delimiter,
            filter_string = "--regex '1:rs.+' --lt '8:0.05'",
            output_filename = "filter_tsv.unzipped.txt",
            header = true
    }
    call TSV.tsv_select as select_unzipped{
        input:
            tsv_input = tsv_input,
            delimiter = delimiter,
            fields="1,2,3,7-10",
            output_filename = "select_tsv.unzipped.txt",
            header = true
    }
    call TSV.tsv_join as join_unzipped{
        input:
            tsv_input = tsv_input,
            tsv_filter_file = input_filter_tsv,
            delimiter = delimiter,
            key_fields = "1",
            data_fields = "1",
            output_filename = "join_tsv.unzipped.txt",
            header = true
    }
    call TSV.tsv_filter as filter_zipped{
        input:
            tsv_input = gzipped_tsv_input,
            delimiter = delimiter,
            filter_string = "--regex '1:rs.+' --lt '8:0.05'",
            output_filename = "filter_tsv.unzipped.txt",
            header = true
    }
    call TSV.tsv_select as select_zipped{
        input:
            tsv_input = gzipped_tsv_input,
            delimiter = delimiter,
            fields="1,2,3,7-10",
            output_filename = "select_tsv.unzipped.txt",
            header = true
    }
    call TSV.tsv_join as join_zipped{
        input:
            tsv_input = gzipped_tsv_input,
            tsv_filter_file = gzipped_input_filter_tsv,
            delimiter = delimiter,
            key_fields = "1",
            data_fields = "1",
            output_filename = "join_tsv.unzipped.txt",
            header = true
    }
    output{
        File append_tsv = tsv_append.tsv_output
        File filter_unzipped_tsv = filter_unzipped.tsv_output
        File select_unzipped_tsv = select_unzipped.tsv_output
        File join_unzipped_tsv = join_unzipped.tsv_output
        File filter_zipped_tsv = filter_zipped.tsv_output
        File select_zipped_tsv = select_zipped.tsv_output
        File join_zipped_tsv = join_zipped.tsv_output
    }
}