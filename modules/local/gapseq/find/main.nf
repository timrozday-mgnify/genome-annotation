process GAPSEQ_FIND {
    tag "${meta.id}"
    label 'process_high'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? projectDir + '/containers/gapseq/gapseq-v2.0.1.sif'
        : 'gapseq:2.0.1'}"

    input:
    tuple val(meta), path(genome), path(db)

    output:
    tuple val(meta), path("*-Reactions.tbl"), emit: reactions
    tuple val(meta), path("*-Pathways.tbl"),  emit: pathways
    tuple val(meta), path("*.faa.gz"),        emit: faa,      optional: true
    tuple val(meta), path("*.gff"),           emit: gff,      optional: true
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: '-p all'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -L ${genome} ${prefix}.fa

    cp -rL ${db} ./gapseq_db_local

    gapseq find \\
        -f ./ \\
        -D ./gapseq_db_local \\
        -K ${task.cpus} \\
        ${args} \\
        ${prefix}.fa

    rm -Rf ./gapseq_db_local

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gapseq: \$(cd /opt/gapseq && git describe --tags 2>/dev/null | sed 's/^v//')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: '-p all'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-all-Reactions.tbl
    touch ${prefix}-all-Pathways.tbl
    touch ${prefix}.faa.gz
    touch ${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gapseq: \$(cd /opt/gapseq && git describe --tags 2>/dev/null | sed 's/^v//')
    END_VERSIONS
    """
}
