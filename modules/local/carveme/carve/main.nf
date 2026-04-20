process CARVEME_CARVE {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/carveme:1.6.6--pyhdfd78af_1'
        : 'quay.io/biocontainers/carveme:1.6.6--pyhdfd78af_1'}"

    input:
    tuple val(meta), path(fasta)
    path mediadb

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.xml"), emit: model
    path "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def solver_arg  = args.contains('--solver') ? '' : '--solver scip'
    def mediadb_arg = (mediadb instanceof List) ? '' : "--mediadb ${mediadb}"
    """
    cp -L ${fasta} ${prefix}.input

    carve \\
        ${solver_arg} \\
        -o ${prefix}.xml \\
        ${mediadb_arg} \\
        ${args} \\
        ${prefix}.input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        carveme: \$(python -c "import carveme; print(carveme.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        carveme: \$(python -c "import carveme; print(carveme.__version__)")
    END_VERSIONS
    """
}
