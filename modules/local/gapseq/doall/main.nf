process GAPSEQ_DOALL {
    tag "${meta.id}"
    label 'process_high'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? projectDir + '/containers/gapseq/gapseq-v2.0.1.sif'
        : 'gapseq:2.0.1'}"

    input:
    // Pass medium as [] to let gapseq auto-predict the growth medium
    tuple val(meta), path(genome), path(db), path(medium)

    output:
    tuple val(meta), path("*-Reactions.tbl"),   emit: reactions
    tuple val(meta), path("*-Pathways.tbl"),    emit: pathways
    tuple val(meta), path("*-Transporter.tbl"), emit: transporters
    tuple val(meta), path("*-draft.RDS"),       emit: draft_rds
    tuple val(meta), path("*-draft.xml"),       emit: draft_xml
    tuple val(meta), path("*-medium.csv"),  optional: true,         emit: medium
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.RDS"),     emit: model_rds
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.xml"),     emit: model_xml
    tuple val(meta), path("*.faa.gz"),          emit: faa,  optional: true
    tuple val(meta), path("*.gff"),             emit: gff,  optional: true
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def medium_cmd  = (medium instanceof List) ? "" : "-m ${medium}"
    """
    cp -L ${genome} ${prefix}.fa

    gapseq doall \\
        -f ./ \\
        -D ${db} \\
        -K ${task.cpus} \\
        ${medium_cmd} \\
        ${args} \\
        ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gapseq: \$(cd /opt/gapseq && git describe --tags 2>/dev/null | sed 's/^v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-all-Reactions.tbl
    touch ${prefix}-all-Pathways.tbl
    touch ${prefix}-Transporter.tbl
    touch ${prefix}-draft.RDS
    touch ${prefix}-draft.xml
    touch ${prefix}-medium.csv
    touch ${prefix}.RDS
    touch ${prefix}.xml
    touch ${prefix}.faa.gz
    touch ${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gapseq: \$(cd /opt/gapseq && git describe --tags 2>/dev/null | sed 's/^v//')
    END_VERSIONS
    """
}
