// Timing pilot: nf-core/proteinannotator (InterProScan with every member database inside
// it), UniFIRE, ESMFold (nf-core/proteinfold 2.0.0's image) and Foldseek on 500 Botryllus
// proteins drawn to match the proteome's length spread. Run through `make run` (see the
// Makefile header) so the head sits in tmux on a compute node and `make watch` shows the
// live task table.
//
//   IPS_DATA ── PROTEINANNOTATOR ── UNIFIRE_LITE   (nested nf-core run, its own task table:
//                                                   `make pa-progress`)
//   ESMFOLD ── FOLDSEEK                            (1 GPU, all 500 proteins in one task;
//                                                   then 500 structures vs AFDB Swiss-Prot)
//
// Every task's time is in the trace (params.pilot_dir/pipeline_info/pilot_trace.txt), which
// bin/scale_timings.py scales to the whole proteome.

nextflow.enable.dsl = 2

// Unpack InterProScan 5.59-91.0's data and press its HMMs, once. storeDir keeps it next to
// the tarball, where ../nf-core-proteinannotator expects it, and skips this task when that
// folder exists, so `make sherlock-untar-ips` there and this task are interchangeable.
// The 5.72 unpack took 1 h 13 min.
process IPS_DATA {
    tag "${params.ips_version}"
    storeDir params.ips_store_dir
    // The unpacked data is tens of GB; storeDir already keeps it, so no copy in results/
    publishDir "${params.pilot_dir}/results/ips_data", enabled: false
    cpus 2
    memory '8 GB'
    time '4h'

    input:
    path tarball
    path md5

    output:
    path "interproscan-${params.ips_version}", emit: dir

    script:
    """
    md5sum -c ${md5}
    tar -pxzf ${tarball}
    cd interproscan-${params.ips_version}
    # 5.59 indexes with initial_setup.py, later versions with setup.py
    if [ -f setup.py ]; then python3 setup.py -f interproscan.properties; else python3 initial_setup.py; fi
    """

    stub:
    """
    mkdir -p interproscan-${params.ips_version}/data
    """
}

// nf-core/proteinannotator at the commit ../nf-core-proteinannotator pins. This task is
// that run's head: it submits the HMM searches, s4pred and InterProScan batches as their
// own Slurm jobs. It launches from a fixed folder so a rerun can -resume it.
process PROTEINANNOTATOR {
    tag 'pilot500'
    cpus 2
    memory '8 GB'
    time '46h'

    input:
    path ips_dir
    path fasta

    output:
    path 'pipeline_info/execution_trace_*.txt', emit: trace
    path 'proteinannotator.log'
    path 'interproscan/*.xml', emit: xml

    script:
    def launch = "${params.pilot_dir}/proteinannotator"
    def extra  = params.pa_extra_config ? "-c ${params.pa_extra_config}" : ''
    """
    ips_data=\$(readlink -f ${ips_dir})/data
    fasta=\$(readlink -f ${fasta})
    task_dir=\$PWD
    mkdir -p ${launch}
    cd ${launch}
    printf 'id,fasta\\npilot500,%s\\n' "\$fasta" > samplesheet.csv
    export JAVA_HOME=${params.java_home} PATH=${params.java_home}/bin:\$PATH
    export NXF_OFFLINE=true NXF_OPTS="-Xms1g -Xmx6g" NXF_SINGULARITY_CACHEDIR=${params.pa_img_dir}
    ${params.nextflow_bin} -log ${launch}/nextflow.log run nf-core/proteinannotator \\
        -r ${params.pipeline_commit} -resume -ansi-log false \\
        -profile singularity -c ${params.pa_config} -c ${projectDir}/pilot.config ${extra} \\
        --input samplesheet.csv --outdir ${launch}/results -work-dir ${launch}/work \\
        --pfam_db ${params.pa_db_dir}/Pfam-A.38.2.hmm.gz \\
        --funfam_db ${params.pa_db_dir}/funfam-hmm3-v4_3_0.lib.gz \\
        --nmpfams_db ${params.pa_db_dir}/nmpfamsdb.hmm.gz \\
        --metagroot_db ${params.pa_db_dir}/metagroot.hmm.gz \\
        --interproscan_db \$ips_data \\
        --interproscan_applications ${params.interproscan_apps} \\
        --interproscan_batch_size ${params.interproscan_batch_size}
    cp -r ${launch}/results/pipeline_info \$task_dir/
    cp ${launch}/nextflow.log \$task_dir/proteinannotator.log
    # The pipeline 'ignore's a failed InterProScan batch and still ends "completed", so
    # count the batches: every one must have written its XML.
    mkdir -p \$task_dir/interproscan
    cp ${launch}/results/functional_annotation/interproscan/*/*.xml \$task_dir/interproscan/
    want=\$(( (\$(grep -c '^>' \$fasta) + ${params.interproscan_batch_size} - 1) / ${params.interproscan_batch_size} ))
    got=\$(ls \$task_dir/interproscan/*.xml | wc -l)
    [ "\$got" -eq "\$want" ] || { echo "InterProScan wrote \$got of \$want batch XMLs" >&2; exit 1; }
    """

    stub:
    """
    mkdir -p pipeline_info interproscan
    touch pipeline_info/execution_trace_stub.txt proteinannotator.log interproscan/pilot500_batch0.xml
    """
}

// UniFIRE's lite image on proteinannotator's InterProScan XML: UniRule, ARBA and PIRSR. The
// FASTA headers carry OX=30301, which UniFIRE reads from the XML. Nextflow runs the image
// with `exec`, so the image's entry script is called by its path.
//
// UniFIRE 2026.3 reads the XML namespace InterProScan uses from 5.66 on; 5.59 writes the
// older one, and UniFIRE then reads every protein without its sequence and stops on each
// system ("this.md5 is null") while still exiting 0. So the batches are merged into one
// file under the newer namespace, which also compiles each rule set once instead of once
// per batch (ARBA's compile took 14 min). Any ERROR line in the log fails the task.
process UNIFIRE_LITE {
    tag 'pilot500'
    container params.unifire_sif
    containerOptions '--writable-tmpfs'
    cpus 4
    memory '32 GB'
    time '6h'

    input:
    path xmls

    output:
    path 'out/predictions*'
    path 'unifire.log'

    script:
    """
    set -o pipefail
    {
        echo '<?xml version="1.0" encoding="UTF-8"?>'
        echo '<protein-matches xmlns="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/schemas" interproscan-version="${params.ips_version}">'
        for xml in ${xmls}; do sed -e '/<protein-matches/d' -e '/<\\/protein-matches>/d' \$xml; done
        echo '</protein-matches>'
    } > pilot500.xml
    mkdir -p out
    /opt/scripts/bin/unifire-workflow.sh -i \$PWD/pilot500.xml -t iprscanxml -o \$PWD/out 2>&1 | tee unifire.log
    if grep -q ' ERROR ' unifire.log; then echo 'UniFIRE logged errors' >&2; exit 1; fi
    """

    stub:
    """
    mkdir -p out
    touch out/predictions_unirule.out unifire.log
    """
}

// nf-core/proteinfold 2.0.0 runs one ESMFold task per protein. Here one GPU task folds all
// 500 (--max-tokens-per-batch 0 folds one protein at a time, as the pipeline does), so the
// GPU is not requested 500 times. The log holds each protein's fold time and the model
// load time, which the scaling adds back once per protein.
process ESMFOLD {
    tag 'pilot500'
    container params.esmfold_sif
    containerOptions '--nv --env TRITON_CACHE_DIR=/tmp/triton_cache --env XDG_CACHE_HOME=/tmp'
    queue 'gpu'
    clusterOptions '-G 1 -C GPU_MEM:32GB|GPU_MEM:48GB|GPU_MEM:80GB'
    cpus 8
    memory '64 GB'
    time '8h'

    input:
    path fasta
    path weights

    output:
    path 'esm-fold.log'
    path 'gpu.txt'
    path 'pdb', emit: pdb

    script:
    """
    set -o pipefail
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader > gpu.txt
    mkdir -p pdb
    esm-fold -i ${fasta} -o pdb -m ${weights} --num-recycles 4 --max-tokens-per-batch 0 \\
        2>&1 | tee esm-fold.log
    """

    stub:
    """
    mkdir -p pdb
    touch esm-fold.log gpu.txt pdb/stub.pdb
    """
}

// Foldseek 9.427df8a, the version nf-core/proteinfold 2.0.0 pins, searching every ESMFold
// structure against AlphaFold DB's Swiss-Prot set (`make stage` downloads it). Proteinfold
// runs this as `foldseek easy-search` per protein; here one task searches all 500, and the
// trace's realtime is the whole search.
process FOLDSEEK {
    tag 'pilot500'
    container params.foldseek_sif
    cpus 16
    memory '64 GB'
    time '8h'

    input:
    path pdb
    path db_dir

    output:
    path 'pilot500_afdb_swissprot.m8'
    path 'foldseek.log'

    script:
    """
    set -o pipefail
    foldseek easy-search ${pdb} ${db_dir}/${params.foldseek_db_name} pilot500_afdb_swissprot.m8 \
        \${L_SCRATCH:-\$PWD}/foldseek-tmp --threads ${task.cpus} \
        --format-output query,target,fident,alnlen,evalue,bits,prob,alntmscore,qtmscore,ttmscore,theader \
        2>&1 | tee foldseek.log
    """

    stub:
    """
    touch pilot500_afdb_swissprot.m8 foldseek.log
    """
}

// proteinannotator removes duplicates by name, keeping one protein per name without an
// error, and ESMFold names each structure file after its protein. So stop before any task
// if a name repeats; `make sample` writes the FASTAs from bin/unique_names.py's output.
def repeated_names(path) {
    file(path, checkIfExists: true).readLines()
        .findAll { it.startsWith('>') }
        .collect { it.substring(1).split(/[\s|]/)[0] }
        .countBy { it }
        .findAll { name, n -> n > 1 }
        .keySet()
}

workflow {
    [params.fasta, params.fasta_ox].each { fa ->
        def dup = repeated_names(fa)
        if (dup) error "${fa}: ${dup.size()} protein name(s) used more than once (${dup.take(3).join(', ')}); run `make sample`"
    }
    ips = IPS_DATA(file(params.ips_tarball, checkIfExists: true),
                   file("${params.ips_tarball}.md5", checkIfExists: true)).dir
    // The FASTA with OX=30301 headers, so UniFIRE can read the taxon from the XML
    PROTEINANNOTATOR(ips, file(params.fasta_ox, checkIfExists: true))
    UNIFIRE_LITE(PROTEINANNOTATOR.out.xml)
    if (!params.skip_fold) {
        ESMFOLD(file(params.fasta, checkIfExists: true), file(params.esmfold_weights, checkIfExists: true))
        FOLDSEEK(ESMFOLD.out.pdb, file(params.foldseek_db_dir, checkIfExists: true))
    }
}
