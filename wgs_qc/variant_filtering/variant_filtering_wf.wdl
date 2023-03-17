task variant_filtering{
    File vcf_in
    String out_prefix

    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 2
    Int max_retries = 3

    command <<<

        # Create IDs for variants, remove singletons variants
        gunzip -c ${vcf_in} |
            perl -lne '
                if (/^#/) {
                    print;
                } else {
                    /;AC=(\d+)/;
                    if ($1 > 1) {
                        /^chr(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)/;
                        $id = join(":", $1, $2, $3, $4);
                        /^(\S+\s+\S+\s+)\S+(\s+.+)/;
                        print $1.$id.$2;
                    }
                }
        ' > ${out_prefix}.non_singletons.vcf

        # Get list of biallelic variants
        perl -lne '
            BEGIN {
                %alleleCounts = ();
                %xref = ();
            }
            if (!/^#/) {
                /^\S+\s+(\S+)\s+(\S+)/;
                if (exists($alleleCounts{$1})) {
                    $alleleCounts{$1}++;
                } else {
                    $alleleCounts{$1} = 1;
                }
                $xref{$1} = $2;
                END {
                    foreach $position (keys(%alleleCounts)) {
                        if ($alleleCounts{$position} == 1) {
                            print $xref{$position};
                        }
                    }
                }
            }
        ' ${out_prefix}.non_singletons.vcf > \
          ${out_prefix}.non_singletons.extract

        # Extract biallelic variants
        perl -lne '
            use warnings;
            use Data::Dumper;
            BEGIN {
                %extract = ();
                open(EXTRACT, "${out_prefix}.non_singletons.extract");
                while (<EXTRACT>) {
                    chomp;
                    $extract{$_} = 1;
                }
                close EXTRACT;
            }
            if (/^#/) {
                print;
            } else {
                /^\S+\s+\S+\s+(\S+)/;
                if (exists($extract{$1})) {
                    print;
                }
            }
        ' ${out_prefix}.non_singletons.vcf > \
          ${out_prefix}.non_singletons.biallelic.vcf

        # Extract only the SVM=PASS variants
        perl -lne '
            if (/^#/) {
                print;
            } elsif (/PASS/) {
                print;
            }
        ' ${out_prefix}.non_singletons.biallelic.vcf > \
          ${out_prefix}.non_singletons.biallelic.svm_pass.vcf

    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File vcf_out = "${out_prefix}.non_singletons.biallelic.svm_pass.vcf"
    }
}


workflow variant_filtering_wf{
    Array[File] vcfs_in
    Array[String] chrs
    Array[String] out_prefixes

    # Call variant filtering task on each chromosome in parallel
    scatter(index in range(length(vcfs_in))){

        call variant_filtering{
            input:
                vcf_in = vcfs_in[index],
                out_prefix = out_prefixes[index],
        }
    }

}



