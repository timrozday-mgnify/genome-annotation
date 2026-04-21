#!/usr/bin/env bash
# Download gapseq reference databases using a local container image.
# Run from anywhere; paths resolve relative to this script's location.
#
# Usage:
#   bash scripts/download_gapseq_db.sh --docker   [--outdir DIR] [-t TAXON]
#   bash scripts/download_gapseq_db.sh --singularity [--outdir DIR] [-t TAXON]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

DOCKER_TAR="${REPO_ROOT}/containers/gapseq/gapseq-v2.0.1.tar"
SINGULARITY_SIF="${REPO_ROOT}/containers/gapseq/gapseq-v2.0.1.sif"
DOCKER_IMAGE="gapseq:2.0.1"

engine=""
outdir="$(pwd)/gapseq_db"
taxon="Bacteria"

usage() {
    cat <<EOF
Usage: bash $(basename "$0") [OPTIONS]

Download gapseq reference databases using a local container image.

Options:
  --docker        Use Docker (loads image from containers/gapseq/gapseq-v2.0.1.tar)
  --singularity   Use Singularity (uses containers/gapseq/gapseq-v2.0.1.sif)
  --outdir DIR    Directory in which to create gapseq_db/ (default: ./gapseq_db)
  -t TAXON        Taxonomy group to download (default: Bacteria; e.g. Archaea)
  -h, --help      Show this help message

After downloading, set in nextflow.config:
  params.databases.gapseq.path = '/absolute/path/to/gapseq_db'
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --docker)      engine="docker"; shift ;;
        --singularity) engine="singularity"; shift ;;
        --outdir)      outdir="$2"; shift 2 ;;
        -t)            taxon="$2"; shift 2 ;;
        -h|--help)     usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [[ -z "$engine" ]]; then
    echo "ERROR: specify --docker or --singularity" >&2
    usage >&2
    exit 1
fi

outdir="$(mkdir -p "$outdir" && cd "$outdir" && pwd)"
db_path="${outdir}/gapseq_db"

echo "=== GapSeq database download ==="
echo "  Engine:  ${engine}"
echo "  Taxon:   ${taxon}"
echo "  Output:  ${db_path}"
echo ""

if [[ "$engine" == "docker" ]]; then
    if [[ ! -f "$DOCKER_TAR" ]]; then
        echo "ERROR: Docker image archive not found: ${DOCKER_TAR}" >&2
        exit 1
    fi
    if [[ -z "$(docker images -q "$DOCKER_IMAGE" 2>/dev/null)" ]]; then
        echo "Loading Docker image from ${DOCKER_TAR} ..."
        docker load -i "$DOCKER_TAR"
        echo ""
    else
        echo "Docker image ${DOCKER_IMAGE} already loaded."
        echo ""
    fi
    docker run --rm \
        --platform linux/amd64 \
        --entrypoint "" \
        -v "${outdir}:/output" \
        "$DOCKER_IMAGE" \
        gapseq update-sequences -t "$taxon" -D /output/gapseq_db

elif [[ "$engine" == "singularity" ]]; then
    if [[ ! -f "$SINGULARITY_SIF" ]]; then
        echo "ERROR: Singularity image not found: ${SINGULARITY_SIF}" >&2
        exit 1
    fi
    singularity exec \
        "$SINGULARITY_SIF" \
        gapseq update-sequences -t "$taxon" -D "${db_path}"
fi

echo ""
echo "=== Done ==="
echo "Database downloaded to: ${db_path}"
echo ""
echo "Add to nextflow.config:"
echo "  params.databases.gapseq.path = '${db_path}'"
