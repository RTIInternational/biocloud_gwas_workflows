import "biocloud_wdl_tools/split_utils/split_utils.wdl" as SPLIT

task array_wc{
    Array[File] inputs
    # Runtime environment
    String docker = "rtibiocloud/pigz:v2.4-8d966cb"
    Int cpu = 4
    Int mem_gb = 8

    command<<<

        # Loop through files in input file and copy/decompress them to output dir
        for input_file in ${sep=" " inputs}; do

            if [[ $input_file == *.gz ]]; then
                pigz -p ${cpu} -d -c $input_file | wc -l
            else
                # Just copy flat files to MultiQC dir
                wc -l $input_file
            fi
        done

    >>>

    output {
        Array[String] wcs = read_lines(stdout())
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }
}


workflow test_split_utils{
    File vcf_info
    File vcf_gz
    File vcf_info_gz
    Int records_per_split = 10000
    Int cpu = 4
    Int mem_gb = 8

    call SPLIT.split_vcf_info as split_info{
        input:
            input_vcf_info = vcf_info,
            records_per_split = records_per_split,
            output_basename = "test_split_vcf",
            cpu = cpu,
            mem_gb = mem_gb
    }


    call SPLIT.split_vcf as split_vcf_gz{
        input:
            input_vcf = vcf_gz,
            records_per_split = records_per_split,
            output_basename = "test_split_vcf",
            cpu = cpu,
            mem_gb = mem_gb
    }

    call SPLIT.split_vcf_info as split_info_gz{
        input:
            input_vcf_info = vcf_info_gz,
            records_per_split = records_per_split,
            output_basename = "test_split_vcf",
            cpu = cpu,
            mem_gb = mem_gb
    }

    call SPLIT.split_vcf_info as split_info_nogzout{
        input:
            input_vcf_info = vcf_info,
            records_per_split = records_per_split,
            output_basename = "test_split_vcf",
            compress_outputs = false,
            cpu = cpu,
            mem_gb = mem_gb
    }

    call SPLIT.split_vcf as split_vcf_gz_nogzout{
        input:
            input_vcf = vcf_gz,
            records_per_split = records_per_split,
            output_basename = "test_split_vcf",
            compress_outputs = false,
            cpu = cpu,
            mem_gb = mem_gb
    }

    call SPLIT.split_vcf_info as split_info_gz_nogzout{
        input:
            input_vcf_info = vcf_info_gz,
            records_per_split = records_per_split,
            output_basename = "test_split_vcf",
            compress_outputs = false,
            cpu = cpu,
            mem_gb = mem_gb
    }

    call array_wc as split_info_wc{
        input:
            inputs = split_info.split_vcf_infos
    }

    call array_wc as split_vcf_gz_wc{
        input:
            inputs = split_vcf_gz.split_vcfs
    }

    call array_wc as split_info_gz_wc{
        input:
            inputs = split_info_gz.split_vcf_infos
    }

    call array_wc as split_info_nogzout_wc{
        input:
            inputs = split_info_nogzout.split_vcf_infos
    }

    call array_wc as split_vcf_gz_nogzout_wc{
        input:
            inputs = split_vcf_gz_nogzout.split_vcfs
    }

    call array_wc as split_info_gz_nogzout_wc{
        input:
            inputs = split_info_gz_nogzout.split_vcf_infos
    }

    output{
        Array[File] split_infos = split_info.split_vcf_infos
        Array[File] split_gz_vcfs = split_vcf_gz.split_vcfs
        Array[File] split_gz_infos = split_info_gz.split_vcf_infos
        Array[File] split_infos_nogzout = split_info_nogzout.split_vcf_infos
        Array[File] split_gz_vcfs_nogzout = split_vcf_gz_nogzout.split_vcfs
        Array[File] split_gz_infos_nogzout = split_info_gz_nogzout.split_vcf_infos

        Array[String] split_infos_wc = split_info_wc.wcs
        Array[String] split_gz_vcfs_wc = split_vcf_gz_wc.wcs
        Array[String] split_gz_infos_wc = split_info_gz_wc.wcs
        Array[String] split_infos_nogzout_wc = split_info_nogzout_wc.wcs
        Array[String] split_gz_vcfs_nogzout_wc = split_vcf_gz_nogzout_wc.wcs
        Array[String] split_gz_infos_nogzout_wc = split_info_gz_nogzout_wc.wcs
    }
}