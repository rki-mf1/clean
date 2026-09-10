// Whether one sample consists of a pair of files.
//
// This is a function and not a derived parameter because it is needed inside
// the process scripts: `addParams()` on include statements was removed with the
// strict syntax of Nextflow >=25.10, and assigning to the params scope from the
// entry workflow is flagged as well.
def lib_pairedness() {
    params.input_type == 'illumina' ? 'paired' : 'single'
}
