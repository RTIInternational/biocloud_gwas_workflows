task make_structure_param_files{
    Int markernames
    Int pop_flag
    Int pop_data
    Int use_pop_info
    Int pfrompopflagonly
    Int label = 1
    Int randomize = 0
    Int burnin = 10000
    Int numreps = 20000
    Int extracols = 0
    Int phased = 0
    Int phaseinfo = 0
    Int noadmix = 0
    Int linkage = 0
    Int locprior = 0

    # Runtime environment
    String docker = "rtibiocloud/structure:v2.3.4-f2d7e82"
    Int cpu = 1
    Int mem_gb = 1

    command <<<
        set -e

        mkdir output
        cp /opt/data/mainparams mainparams
        cp /opt/data/extraparams extraparams

        #### UPDATE MAINPARAMS
        # Update mainparams to indicate popflag in use
        echo "Updating popflag..."
        sed -i 's/#define POPFLAG\s*[0|1]/#define POPFLAG   ${pop_flag}/g' mainparams

        echo "Updating popdata..."
        sed -i 's/#define POPDATA\s*[0|1]/#define POPDATA   ${pop_data}/g' mainparams

        # Update mainparams to indicate whether header row
        echo "Updating markernames..."
        sed -i 's/#define MARKERNAMES\s*[0|1]/#define MARKERNAMES      ${markernames}/g' mainparams

        # Update mainparams to indicate whether header row
        echo "Updating label..."
        sed -i 's/#define LABEL\s*[0|1]/#define LABEL      ${label}/g' mainparams

        # Update burn-in
        echo "Updating burnin..."
        sed -i 's/#define BURNIN\s*[0-9]*/#define BURNIN    ${burnin}/g' mainparams

        # Update numreps
        echo "Updating numreps..."
        sed -i 's/#define NUMREPS\s*[0-9]*/#define NUMREPS   ${numreps}/g' mainparams

        # Update phased
        echo "Updating phased..."
        sed -i 's/#define PHASED\s*[0|1]/#define PHASED           ${phased}/g' mainparams

        # Update extracols
        echo "Updating extracols..."
        sed -i 's/#define EXTRACOLS\s*[0|1]/#define EXTRACOLS ${extracols}/g' mainparams

        # Update phaseinfo
        echo "Updating phaseinfo..."
        sed -i 's/#define PHASEINFO\s*[0|1]/#define PHASEINFO        ${phaseinfo}/g' mainparams

        #### UPDATE EXTRAPARAMS
        # Update whether to use pop_info
        echo "Updating usepopinfo..."
        sed -i 's/#define USEPOPINFO\s*[0|1]/#define USEPOPINFO      ${use_pop_info}/g' extraparams

        # Update extraparams to not randomize seed
        echo "Updating randomize..."
        sed -i 's/#define RANDOMIZE\s*[0|1]/#define RANDOMIZE      ${randomize}/g' extraparams

        # Update extraparams with whether to use admixture model
        echo "Updating noadmix..."
        sed -i 's/#define NOADMIX\s*[0|1]/#define NOADMIX      ${noadmix}/g' extraparams

        # Update extraparams with whether to use linkage model
        echo "Updating linkage..."
        sed -i 's/#define LINKAGE\s*[0|1]/#define LINKAGE      ${linkage}/g' extraparams

        # Update extraparams with whether to use location as a prior
        echo "Updating locprior..."
        sed -i 's/#define LOCPRIOR\s*[0|1]/#define LOCPRIOR      ${locprior}/g' extraparams

        echo "Updating pfrompopflagonly..."
        sed -i 's/#define PFROMPOPFLAGONLY\s*[0|1]/#define PFROMPOPFLAGONLY      ${pfrompopflagonly}/g' extraparams
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File mainparams_out = "mainparams"
        File extraparams_out = "extraparams"
    }
}

task structure{
    File input_file
    File mainparams
    File extraparams
    File? stratparams
    String output_basename

    Int k
    Int numloci

    Int? seed
    Int default_seed = 1523031945
    Int actual_seed = select_first([seed, default_seed])

    # Runtime environment
    String docker = "rtibiocloud/structure:v2.3.4-f2d7e82"
    Int cpu = 8
    Int mem_gb = 16

    command <<<
        set -e

        mkdir structure_output

        # Count number of individuals
        num_inds=$(wc -l ${input_file} | perl -lane 'print $F[0]/2;')

        # Create symlink for input bc STRUCTURE is stoopid and will throw a buffer error for long cromwell paths
        ln -s ${input_file} structure.in

        structure -K ${k} \
            -m ${mainparams} \
            ${'-e ' + extraparams} \
            ${'-s ' + stratparams} \
            -L ${numloci} \
            -N $num_inds \
            -D ${default_seed} \
            -i structure.in \
            -o structure_output/${output_basename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File structure_out = "structure_output/${output_basename}_f"
    }
}

task ped2structure{
    # Utility for converting ped file to STRUCTURE-formatted input
    File ped_in
    String output_filename

    # Optional file of reference samples that should have pop info and pop_flag set
    File? ref_samples
    Int? ref_pop_col = 2
    String? ref_delim = "space"

    # Runtime environment
    String docker = "rtibiocloud/ped2structure:v1.0-c3278c6"
    Int cpu = 1
    Int mem_gb = 2

    command<<<
        python /opt/ped2structure.py --ped ${ped_in} \
            ${'--ref-samples ' + ref_samples} \
            ${'--ref-pop-col ' + ref_pop_col} \
            ${'--ref-delim ' + ref_delim} \
            --output ${output_filename} \
            -vvv
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File structure_input = "${output_filename}"
    }
}

task parse_structure_output{
    # Utility for converting ped file to STRUCTURE-formatted input
    File structure_output
    String output_filename
    String delim = " "

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command<<<
        set -e

        # Parse ancestry proportions from structure output
        perl -ne '
          if (/%Miss/) {$in=1;}
          if ($in==1 && !/Label/ && !/^\s+$/) {
              # Remove leading spaces
              s/^\s+//g;

              # Split on whitespace
              @F=split /\s+/;

              if ($F[6] =~ m/\|/){
                # Case: Handle ref samples with pipes in admix proportions
                # Count total number of pops
                @pop_count = (join " ", @F) =~ /\|/g;
                $num_pops = @pop_count;
                # Set the corresponding admix proportion for sample to 1.000, rest to 0.000
                $pop = $F[3];
                @pop_freqs = ("0.000") x $num_pops;
                $pop_freqs[($pop-1)] = "1.000";
                $data = join "${delim}", @pop_freqs;
              }
              else{
                # Case: normal samples with actual admix proportions
                # Grab all admixture proportions (index 5 to the end)
                $data = join "${delim}", splice @F, 5;
              }

              # Remove parentheses from missing data pct
              $F[2] =~ y/()//d;

              # Add STRUCTURE pop label
              $data = $F[0]."${delim}".$F[1]."${delim}".$F[2]." ".$F[3]."${delim}".$data;
              print $data."\n";
          }
          s/\s+//g;
          if ($_ eq "") { $in=0; }' ${structure_output} > ${output_filename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File structure_out = "${output_filename}"
    }
}
