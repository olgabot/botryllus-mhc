// Timing pilot: nf-core/proteinannotator, ESMFold (nf-core/proteinfold 2.0.0's image),
// InterProScan with every member database, and UniFIRE on 500 Botryllus proteins drawn to
// match the proteome's length spread. Run through `make run` (see the Makefile header) so
// the head sits in tmux on a compute node and `make watch` shows the live task table.
//
//   IPS_DATA ──┬── PROTEINANNOTATOR            (nested nf-core run, its own task table:
//              │                                 `make pa-progress`)
//              └── INTERPROSCAN_FULL ── UNIFIRE_LITE
//   ESMFOLD                                     (1 GPU, all 500 proteins in one task)
//
// Every task's time is in the trace (params.pilot_dir/pipeline_info/pilot_trace.txt), which
// bin/scale_timings.py scales to the whole proteome.

nextflow.enable.dsl = 2

// Unpack InterProScan 5.59-91.0's data and press its HMMs, once. storeDir keeps it next to
// the tarball, where ../nf-core-proteinannotator expects it. The 5.72 unpack took 1 h 13 min.
process IPS_DATA {
    tag "${params.ips_version}"
    storeDir params.ips_store_dir
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
    python3 setup.py -f interproscan.properties || python3 setup.py interproscan.properties
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
    time '2d'

    input:
    path ips_dir
    path fasta

    output:
    path 'pipeline_info/execution_trace_*.txt', emit: trace
    path 'proteinannotator.log'

    script:
    def launch = "${params.pilot_dir}/proteinannotator"
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
        -profile singularity -c ${params.pa_config} -c ${projectDir}/pilot.config \\
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
    """

    stub:
    """
    mkdir -p pipeline_info
    touch pipeline_info/execution_trace_stub.txt proteinannotator.log
    """
}

// Every member database InterProScan 5.59 ships (no --applications), set up the way the
// proteinannotator module does it: copy interproscan.properties, point bin.directory at
// the image. The FASTA carries OX=30301 because UniFIRE reads the taxon from this XML.
process INTERPROSCAN_FULL {
    tag 'pilot500'
    container params.ips_image
    cpus 16
    memory '32 GB'
    time '24h'

    input:
    path ips_dir
    path fasta

    output:
    path 'pilot500_full.xml', emit: xml
    path 'pilot500_full.tsv'

    script:
    """
    ln -sfn \$(readlink -f ${ips_dir})/data data
    INTERPROSCAN_DIR="\$(dirname "\$(dirname "\$(which interproscan.sh)")")"
    cp "\$(find "\$INTERPROSCAN_DIR/share" -name interproscan.properties)" .
    sed -i "/^bin\\.directory=/ s|.*|bin.directory=\$INTERPROSCAN_DIR/bin|" interproscan.properties
    export INTERPROSCAN_CONF=interproscan.properties
    tmp=\${L_SCRATCH:-\$PWD}/ips-tmp
    interproscan.sh --cpu ${task.cpus} --disable-precalc --input ${fasta} \\
        --formats XML,TSV --output-file-base pilot500_full --tempdir \$tmp
    """

    stub:
    """
    touch pilot500_full.xml pilot500_full.tsv
    """
}

// UniFIRE's lite image on that XML: UniRule, ARBA and PIRSR. Nextflow runs the image with
// `exec`, so the image's entry script is called by its path.
process UNIFIRE_LITE {
    tag 'pilot500'
    container params.unifire_sif
    containerOptions '--writable-tmpfs'
    cpus 4
    memory '32 GB'
    time '6h'

    input:
    path xml

    output:
    path 'out/predictions*'
    path 'unifire.log'

    script:
    """
    set -o pipefail
    mkdir -p out
    /opt/scripts/bin/unifire-workflow.sh -i \$PWD/${xml} -t iprscanxml -o \$PWD/out 2>&1 | tee unifire.log
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
    touch esm-fold.log gpu.txt
    """
}

workflow {
    ips = IPS_DATA(file(params.ips_tarball, checkIfExists: true),
                   file("${params.ips_tarball}.md5", checkIfExists: true)).dir
    PROTEINANNOTATOR(ips, file(params.fasta, checkIfExists: true))
    INTERPROSCAN_FULL(ips, file(params.fasta_ox, checkIfExists: true))
    UNIFIRE_LITE(INTERPROSCAN_FULL.out.xml)
    ESMFOLD(file(params.fasta, checkIfExists: true), file(params.esmfold_weights, checkIfExists: true))
}
