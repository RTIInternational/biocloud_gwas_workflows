workflow genome_liftover {
    String chainfile =  "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"

    String docker = "python:3.8" # python3.8 has wget tool

    call download_chainfile {
        input:
            liftover_file = liftover_file,
            docker = docker
    }
}

task download_chainfile {
    String chainfile
    String outfile = basename(chainfile, ".gz") 
    #String remove_path = sub(liftover_file, ".+\\/", "") 

    String docker
    Int cpu = 1
    Int mem = 4

    command {
        wget ${chainfile}
        gunzip ${outfile}".gz"
    }

    output {
        File unzipped_output = "${outfile}"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem} GB"
    }

    parameter_meta {
        liftover_file: "Liftover file. Navigate to the webpage http://hgdownload.soe.ucsc.edu/downloads.html and from there click 'LiftOver files' under the  Human genome section. If my starting build is hg38, I would look under the Dec. 2013 (GRCh38/hg38) subsection. Clicking 'LiftOver files' will take you to a page that contains the UCSC Chain Files. Select the liftover file that converts to the build you desire."
        docker: "Docker image."
        cpu: "Number of CPUs for the image."
        mem: "Amount of RAM in GB for the image."
    }

    meta {
        description: "Download and decompress the liftover file."
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}
