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

####################################################################################################

task create_bedfile {

    File input_sumstats
    Int chromosome_col
    Int position_col

    String docker
    Int cpu = 1
    Int mem = 2

    command 
    <<<

        python - <<EOF
        import gzip

        infile = "${input_sumstats}"
        outfile = "sumstats.bed"


        with gzip.open(infile, 'rt') as inF, open(outfile, 'w') as outF:
            head = inF.readline()
            line = inF.readline()

            print(head)
            print(line)
            while line:
                sl = line.split()
                chrom = "chr{}".format(sl[${chromosome_col}])
                pos = sl[${position_col}]
                pos_plus = str(int(pos) + 1)
                outline = "{}\t{}\t{}\n".format(chrom, pos, pos_plus)
                outF.write(outline)
                line = inF.readline()
        EOF
    >>>

    output {
        File sumstats_bed = "sumstats.bed"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem} GB"
    }

    parameter_meta {
        docker: "Docker image."
        cpu: "Number of CPUs for the image."
        mem: "Amount of RAM in GB for the image."
    }

    meta {
        description: "Create BED formatted file from the input summary stats. See BED format: http://genome.ucsc.edu/FAQ/FAQformat.html#format1"
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}

####################################################################################################

task perform_liftover {

    File input_bedfile
    String output_bedfile = "new_build.bed"
    String unmapped_bedfile = "unmapped.bed"
    File chainfile

    String docker
    Int cpu = 1
    Int mem = 2

    command 
    <<<
    # liftOver - Move annotations from one assembly to another
    # usage:
    #   liftOver oldFile map.chain newFile unMapped
    /opt/liftOver \
        ${input_bedfile}  ${chainfile} \
        ${output_bedfile}  ${unmapped_bedfile}

    >>>

    output {
        File output_bed = "${output_bedfile}"
        File unmapped_bed = "${unmapped_bedfile}"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem} GB"
    }

    parameter_meta {

        input_bedfile: ""
        output_bedfile: ""
        unmapped_bedfile: ""
        chainfile: ""
        docker: "Docker image."
        cpu: "Number of CPUs for the image."
        mem: "Amount of RAM in GB for the image."
    }

    meta {
        description: "Perform liftover on the BED formatted sumstats file. See BED format: http://genome.ucsc.edu/FAQ/FAQformat.html#format1"
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}



####################################################################################################
task  final_sumstats {

    File original_sumstats
    String new_build
    File new_bed
    File unmapped_bed
    String output_name
    String unmapped_name = "unmapped_sumstats.txt"
    Int position_col
    Int chr_col

    String docker
    Int cpu = 1
    Int mem = 2

    command 
    <<<

    python - <<EOF
    import gzip

    original_sumstats = "${original_sumstats}"
    new_bed = "${new_bed}"
    unmapped = "${unmapped_bed}"
    outfile = "${output_name}"
    out_unmapped = "${unmapped_name}"
    position_col = ${position_col}
    chr_col = ${chr_col}
    new_build = ${new_build}

    with gzip.open(original_sumstats, "rt") as sumSTATS, \
        open(new_bed) as newBED, \
        open(unmapped) as unMAPPED, \
        open(out_unmapped, "w") as outUNMAPPED, \
        open(outfile, "w") as outF:

        line = unMAPPED.readline()
        unmapped_set = set() # keep track of all of the unmapped positions
        while line:
            if line[0] != "#":
                sl = line.split()
                chrom = sl[0][3:]
                pos = sl[1]
                chr_pos = "{}:{}".format(chrom, pos)
                unmapped_set.add(chr_pos) # add the original chr and position that did not mapped to new genome build

            line = unMAPPED.readline()
        
        header = sumSTATS.readline().split()
        print(header)
        outUNMAPPED.write(" ".join(header) + "\n")
        header[position_col] = "hg{}_POS".format(new_build)
        header = " ".join(header) + "\n" # change to space separated
        outF.write(header)
    
        stats_line = sumSTATS.readline()
        newbed_line = newBED.readline()
        while stats_line:
            sl = stats_line.split()
            position = sl[position_col]
            chromosome = sl[chr_col]
            snp_id = "{}:{}".format(chromosome, position)
            if snp_id not in unmapped_set: # if it's not in the unmapped set, update the position and print the line
                newbed_split = newbed_line.split()
                sl[position_col] = newbed_split[1]
                outline = " ".join(sl) + "\n"
                outF.write(outline)
                newbed_line = newBED.readline()
            else: # else if it is in the unmapped set, add it to the unmapped output
                unmapped_line = " ".join(sl) + "\n"
                outUNMAPPED.write(unmapped_line)
    
            stats_line = sumSTATS.readline()
    
    EOF
    >>>

    output {
        File output_sumstats = "${output_name}"
        File unmapped_sumstats = "${unmapped_name}"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem} GB"
    }

    parameter_meta {
        original_sumstats: "Original set of summary statistics."
        new_bed: "Variants mapped to new build: BED file."
        unmapped_bed: "Variants not mapped to new build: BED file."
        output_name: "Name of final sumstats file with build liftover."
        position_col: "0-based column of position."
        chr_col: "0-based column of chromosome."
        docker: "Docker image."
        cpu: "Number of CPUs for the image."
        mem: "Amount of RAM in GB for the image."
    }

    meta {
        description: ""
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}
