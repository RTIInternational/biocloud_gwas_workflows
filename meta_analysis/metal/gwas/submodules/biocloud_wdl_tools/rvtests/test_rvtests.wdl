import "biocloud_wdl_tools/rvtests/rvtests.wdl" as RVTESTS

workflow test_rvtests{
    call RVTESTS.rvtests
    output{
        Array[File] outputs = rvtests.outputs
    }
}