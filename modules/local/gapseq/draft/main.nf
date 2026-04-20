process GAPSEQ_DRAFT {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? projectDir + '/containers/gapseq/gapseq-v2.0.1.sif'
        : 'gapseq:2.0.1'}"

    input:
    tuple val(meta), path(reactions), path(transporters), path(pathways)

    output:
    tuple val(meta), path("*-draft.RDS"), emit: draft_rds
    tuple val(meta), path("*-draft.xml"), emit: draft_xml
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gapseq draft \\
        -r ${reactions} \\
        -t ${transporters} \\
        -p ${pathways} \\
        -n ${prefix} \\
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
    touch ${prefix}-draft.RDS
    touch ${prefix}-draft.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gapseq: \$(cd /opt/gapseq && git describe --tags 2>/dev/null | sed 's/^v//')
    END_VERSIONS
    """
}
