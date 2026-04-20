process GAPSEQ_FILL {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? projectDir + '/containers/gapseq/gapseq-v2.0.1.sif'
        : 'gapseq:2.0.1'}"

    input:
    tuple val(meta), path(draft_model, stageAs: "input_model/*"), path(medium, stageAs: "input_medium/*")

    output:
    tuple val(meta), path("*.RDS"), emit: model_rds
    tuple val(meta), path("*.xml"), emit: model_xml
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gapseq fill \\
        -m ${draft_model} \\
        -n ${medium} \\
        -f ./ \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gapseq: \$(cd /opt/gapseq && git describe --tags 2>/dev/null | sed 's/^v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.RDS
    touch ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gapseq: \$(cd /opt/gapseq && git describe --tags 2>/dev/null | sed 's/^v//')
    END_VERSIONS
    """
}
