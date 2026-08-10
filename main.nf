#!/usr/bin/env nextflow

/*
 * digestome: licence-clean functional profiling of anaerobic digestion microbiomes.
 *
 * This workflow covers the ANALYSIS path only: proteomes in, per-genome scores and
 * one community profile out. Database staging and panel construction stay as the
 * scripts in scripts/, because they are one-time, need internet access, and gain
 * nothing from being a DAG.
 *
 * The processes call the same command-line tools a user would run by hand. That is
 * deliberate: the tools remain usable standalone with nothing but HMMER and a
 * Python interpreter, and this workflow is a wrapper for scale and portability
 * rather than a replacement for them.
 */

nextflow.enable.dsl = 2

def helpMessage() {
    log.info """
    digestome ${workflow.manifest.version}

    Usage:
      nextflow run . --proteomes 'path/to/*.faa' --db /path/to/databases [options]

    Required:
      --proteomes     Glob of protein FASTA files, one per genome/MAG. Quote it.
                      Alternatively --input a CSV with columns: sample,faa
      --db            Database root containing ad_panel/ (from scripts/build_ad_panel.sbatch)

    Optional:
      --gtdbtk        GTDB-Tk summary.tsv. Strongly recommended: without it the
                      acetoclastic call falls back to gene evidence, which over-calls.
      --secretion     Extracellular-targeting modules (default: bundled)
      --panel         Panel TSV (default: bundled)
      --sample_name   Name for the community profile
      --outdir        Results directory (default: results)

    Profiles:
      -profile slurm | conda | singularity | docker | test
    """.stripIndent()
}

process HMMSEARCH {
    tag   "$sample"
    label 'process_low'

    input:
    tuple val(sample), path(faa)
    path  hmm

    output:
    tuple val(sample), path("${sample}.tblout")

    script:
    // --cut_nc applies each model's own curated cutoff. Never substitute a global
    // -T: curated cutoffs across this panel span roughly 20 to 1400 bits.
    """
    hmmsearch --cpu ${task.cpus} --cut_nc --tblout ${sample}.tblout ${hmm} ${faa} > /dev/null
    """

    stub:
    "touch ${sample}.tblout"
}

process SCORE {
    tag   "$sample"
    label 'process_single'

    input:
    tuple val(sample), path(tblout)
    path  map
    path  panel
    path  secretion
    path  gtdbtk

    output:
    path "${sample}.modules.tsv", emit: modules
    path "${sample}.summary.txt", emit: summary

    script:
    def sec = secretion.name != 'NO_FILE' ? "--secretion ${secretion}" : ''
    def tax = gtdbtk.name    != 'NO_FILE' ? "--gtdbtk ${gtdbtk}"      : ''
    """
    python3 ${projectDir}/panel_scored.py \\
        --tblout ${tblout} --map ${map} --panel ${panel} \\
        --name ${sample} --out ${sample}.modules.tsv \\
        ${sec} ${tax} > ${sample}.summary.txt
    """

    stub:
    """
    printf 'genome\\tmodule\\tbranch\\tcore_expected\\tcore_detectable\\tcore_found\\tcore_secreted\\tpct_complete\\tstatus\\n' > ${sample}.modules.tsv
    printf '# ${sample} stub\\n' > ${sample}.summary.txt
    """
}

process AGGREGATE {
    label 'process_single'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path scored, stageAs: 'scored/*'

    output:
    path 'profile.json'
    path 'profile.txt'

    script:
    """
    python3 ${projectDir}/aggregate_community.py \\
        --dir scored \\
        --sample '${params.sample_name}' \\
        --description '${params.description}' \\
        --out-json profile.json --out-txt profile.txt
    """

    stub:
    """
    echo '{}' > profile.json
    touch profile.txt
    """
}

workflow {

    // The 26.x strict parser rejects top-level statements, so the help check lives
    // here rather than beside the function.
    if (params.help) {
        helpMessage()
        return
    }

    // ---- inputs -------------------------------------------------------------
    if (!params.proteomes && !params.input) {
        error "Provide --proteomes '<glob>' or --input <samplesheet.csv>. See --help."
    }

    ch_faa = params.input
        ? Channel.fromPath(params.input, checkIfExists: true)
              .splitCsv(header: true)
              .map { row ->
                  if (!row.sample || !row.faa) error "Samplesheet needs columns: sample,faa"
                  tuple(row.sample, file(row.faa, checkIfExists: true))
              }
        : Channel.fromPath(params.proteomes, checkIfExists: true)
              .map { f -> tuple(f.simpleName, f) }

    // Fail early and clearly rather than deep inside hmmsearch.
    def dbdir = params.db ? file("${params.db}/ad_panel") : null
    if (!dbdir || !dbdir.exists()) {
        error "Panel database not found at \${params.db}/ad_panel. Build it first:\n" +
              "  DB=<root> bash scripts/prefetch_ncbifam.sh && bash scripts/prefetch_pfam.sh\n" +
              "  DB=<root> sbatch scripts/build_ad_panel.sbatch"
    }

    ch_hmm   = Channel.value(file("${params.db}/ad_panel/ad_panel.hmm",     checkIfExists: true))
    ch_map   = Channel.value(file("${params.db}/ad_panel/ad_panel_map.tsv", checkIfExists: true))
    ch_panel = Channel.value(file(params.panel,     checkIfExists: true))
    ch_sec   = Channel.value(params.secretion ? file(params.secretion, checkIfExists: true)
                                              : file("${projectDir}/assets/NO_FILE"))
    ch_tax   = Channel.value(params.gtdbtk    ? file(params.gtdbtk,    checkIfExists: true)
                                              : file("${projectDir}/assets/NO_FILE"))

    if (!params.gtdbtk) {
        log.warn "No --gtdbtk supplied: the acetoclastic call will fall back to gene " +
                 "evidence only, which over-calls. ACDS/CODH is reversible and " +
                 "autotrophic hydrogenotrophs carry it for carbon fixation."
    }

    // ---- analysis -----------------------------------------------------------
    HMMSEARCH(ch_faa, ch_hmm)
    SCORE(HMMSEARCH.out, ch_map, ch_panel, ch_sec, ch_tax)

    // One aggregate over every scored genome. The scorer writes two files per
    // genome and the aggregator reads both, so mix them into a single collection.
    AGGREGATE(SCORE.out.modules.mix(SCORE.out.summary).collect())
}
