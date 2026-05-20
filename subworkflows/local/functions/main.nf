/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Remove empty fasta files from a channel
//
def rmEmptyFastAs(ch_fastas) {
    // if not stub run, remove empty fasta files
    if (!workflow.stubRun) {
        def ch_nonempty_fastas = ch_fastas
            .filter { _meta, fasta ->
                try {
                    file(fasta).countFasta( limit: 1 ) > 0
                } catch (java.util.zip.ZipException _e) {
                    log.debug "[rmEmptyFastAs]: ${fasta} is not in GZIP format, this is likely because it was cleaned with --remove_intermediate_files"
                    true
                } catch (_EOFException) {
                    log.debug "[rmEmptyFastAs]: ${fasta} has an EOFException, this is likely an empty gzipped file."
                }
            }
        return ch_nonempty_fastas
    }
    // if stub run, return the input channel
    else {
        return ch_fastas
    }
}

//
// Remove empty tsv files from a channel
//
def rmEmptyTsvs(ch_tsvs) {
    // if not stub run, remove empty tsv files
    if (!workflow.stubRun) {
        def ch_nonempty_tsvs = ch_tsvs
            .filter { _meta, tsv ->
                try {
                    file(tsv).countLines( limit: 2 ) > 1
                } catch (java.util.zip.ZipException _e) {
                    log.debug "[rmEmptyTsvss]: ${tsv} is not in GZIP format, this is likely because it was cleaned with --remove_intermediate_files"
                    true
                } catch (_EOFException) {
                    log.debug "[rmEmptyTsvss]: ${tsv} has an EOFException, this is likely an empty gzipped file."
                }
            }
        return ch_nonempty_tsvs
    }
    // if stub run, return the input channel
    else {
        return ch_tsvs
    }
}

//
// Count the number of sequences in a channel of fasta files
//
def countFastAs(ch_fastas) {
    // if not stub run, count the number of sequences in the fasta files
    if (!workflow.stubRun) {
        def fasta_counts = 0
        ch_fastas
            .map { _meta, fasta ->
                try {
                    fasta_counts += file(fasta).countFasta( limit: 2 )
                } catch (java.util.zip.ZipException _e) {
                    log.debug "[rmEmptyFastAs]: ${fasta} is not in GZIP format, this is likely because it was cleaned with --remove_intermediate_files"
                    true
                } catch (_EOFException) {
                    log.debug "[rmEmptyFastAs]: ${fasta} has an EOFException, this is likely an empty gzipped file."
                }
            }
        return fasta_counts
    }
    // if stub run, return 0
    else {
        return 0
    }
}