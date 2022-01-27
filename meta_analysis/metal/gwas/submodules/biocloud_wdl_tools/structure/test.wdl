import "biocloud_wdl_tools/structure/structure.wdl" as STRUCT

workflow test_struct_params{

    String derp

    call STRUCT.make_structure_param_files{
        input:
            markernames = 1,
            pop_flag = 2,
            use_pop_info = 3,
            pop_data = 4,
            burnin = 5,
            numreps = 6,
            label = 7,
            randomize = 8,
            extracols = 9,
            phased = 10,
            phaseinfo = 11,
            noadmix = 12,
            linkage = 13,
            locprior = 14
    }

    output{
        File mainparams = make_structure_param_files.mainparams_out
        File extraparams = make_structure_param_files.extraparams_out
    }
}