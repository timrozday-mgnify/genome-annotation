#!/usr/bin/env bash
# One-time script to create mini test databases from full databases in ~/Downloads/.
# Run from the pipeline root directory:
#   bash tests/scripts/make_mini_databases.sh
set -eu
# pipefail is intentionally not set — gunzip exits with SIGPIPE when awk finishes early,
# which would fail the whole pipeline. We check outputs explicitly instead.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TESTDATA="${REPO_ROOT}/tests/testdata"

PFAM_SRC="${HOME}/Downloads/Pfam-A.hmm.gz"
DBCAN_SRC="${HOME}/Downloads/dbCAN-HMMdb-V14.txt"
CAZY_FA_SRC="${HOME}/Downloads/CAZyDB.07242025.fa"
KOLIST_SRC="${HOME}/Downloads/ko_list.gz"
PROFILES_SRC="${HOME}/Downloads/profiles.tar.gz"

echo "=== Mini Pfam ==="
mkdir -p "${TESTDATA}/mini_pfam"
gunzip -c "${PFAM_SRC}" 2>/dev/null \
  | awk 'BEGIN{n=0} /^HMMER3/{n++} n<=100{print} n>100{exit}' \
  > "${TESTDATA}/mini_pfam/mini_pfam.hmm" || true
[ -s "${TESTDATA}/mini_pfam/mini_pfam.hmm" ] || { echo "ERROR: mini_pfam.hmm is empty"; exit 1; }

docker run --rm --platform linux/amd64 --entrypoint /bin/bash \
  -v "${TESTDATA}/mini_pfam:/db" \
  quay.io/biocontainers/hmmer:3.4--hb6cb901_4 \
  -c "hmmpress /db/mini_pfam.hmm"

echo "=== Mini dbCAN ==="
mkdir -p "${TESTDATA}/mini_dbcan"
DBCAN_BASE="https://pro.unl.edu/dbCAN2/download/run_dbCAN_database_total"

# Download support files (bcb.unl.edu redirects to HTML; pro.unl.edu has the actual files)
echo "  Downloading and trimming dbCAN-sub.hmm (first 50 subfamilies)..."
curl -sSL "${DBCAN_BASE}/dbCAN_sub.hmm" \
  | awk 'BEGIN{n=0} /^HMMER3/{n++} n<=50{print} n>50{exit}' \
  > "${TESTDATA}/mini_dbcan/dbCAN-sub.hmm" || true
[ -s "${TESTDATA}/mini_dbcan/dbCAN-sub.hmm" ] || { echo "ERROR: dbCAN-sub.hmm is empty"; exit 1; }
echo "  Downloading fam-substrate-mapping.tsv..."
curl -sSL "${DBCAN_BASE}/fam-substrate-mapping.tsv" -o "${TESTDATA}/mini_dbcan/fam-substrate-mapping.tsv"

# Build mini dbCAN.hmm: first 50 of 875 CAZyme families
awk 'BEGIN{n=0} /^HMMER3/{n++} n<=50{print} n>50{exit}' \
  "${DBCAN_SRC}" \
  > "${TESTDATA}/mini_dbcan/dbCAN.hmm" || true
[ -s "${TESTDATA}/mini_dbcan/dbCAN.hmm" ] || { echo "ERROR: dbCAN.hmm is empty"; exit 1; }

# Build mini CAZy.dmnd from first 500 sequences of local FASTA
awk '/^>/{n++} n>500{exit} {print}' "${CAZY_FA_SRC}" > "${TESTDATA}/mini_dbcan/mini_cazy.fa"
[ -s "${TESTDATA}/mini_dbcan/mini_cazy.fa" ] || { echo "ERROR: mini_cazy.fa is empty"; exit 1; }

docker run --rm --platform linux/amd64 --entrypoint /bin/bash \
  -v "${TESTDATA}/mini_dbcan:/db" \
  quay.io/biocontainers/diamond:2.1.9--h43eeafb_0 \
  -c "diamond makedb --in /db/mini_cazy.fa --db /db/CAZy"

# hmmpress dbCAN.hmm and dbCAN-sub.hmm
docker run --rm --platform linux/amd64 --entrypoint /bin/bash \
  -v "${TESTDATA}/mini_dbcan:/db" \
  quay.io/biocontainers/hmmer:3.4--hb6cb901_4 \
  -c "hmmpress /db/dbCAN.hmm && hmmpress /db/dbCAN-sub.hmm"

echo "=== Mini KOfam ==="
mkdir -p "${TESTDATA}/mini_kofam/chunk_1/profiles"
mkdir -p "${TESTDATA}/mini_kofam/chunk_2/profiles"

# Write ko_list files (header + 50 entries each)
gunzip -c "${KOLIST_SRC}" | head -1 \
  > "${TESTDATA}/mini_kofam/chunk_1/ko_list"
gunzip -c "${KOLIST_SRC}" | grep -v "^knum" | head -50 \
  >> "${TESTDATA}/mini_kofam/chunk_1/ko_list"

gunzip -c "${KOLIST_SRC}" | head -1 \
  > "${TESTDATA}/mini_kofam/chunk_2/ko_list"
gunzip -c "${KOLIST_SRC}" | grep -v "^knum" | sed -n '51,100p' \
  >> "${TESTDATA}/mini_kofam/chunk_2/ko_list"

# Extract matching HMM profiles from profiles.tar.gz
echo "  Extracting chunk_1 profiles..."
while IFS=$'\t' read -r ko _rest; do
  tar -xzf "${PROFILES_SRC}" \
    -C "${TESTDATA}/mini_kofam/chunk_1/profiles" \
    --strip-components=1 \
    "profiles/${ko}.hmm" 2>/dev/null || true
done < <(tail -n +2 "${TESTDATA}/mini_kofam/chunk_1/ko_list")

echo "  Extracting chunk_2 profiles..."
while IFS=$'\t' read -r ko _rest; do
  tar -xzf "${PROFILES_SRC}" \
    -C "${TESTDATA}/mini_kofam/chunk_2/profiles" \
    --strip-components=1 \
    "profiles/${ko}.hmm" 2>/dev/null || true
done < <(tail -n +2 "${TESTDATA}/mini_kofam/chunk_2/ko_list")

echo "=== Samplesheet ==="
cat > "${TESTDATA}/samplesheet.csv" <<CSV
sample,fastx
mini_genome,${TESTDATA}/mini_genome.fa
CSV

echo ""
echo "Done. Mini database sizes:"
du -sh "${TESTDATA}/mini_pfam" "${TESTDATA}/mini_dbcan" "${TESTDATA}/mini_kofam"
echo ""
echo "Profile counts:"
echo "  chunk_1: $(ls ${TESTDATA}/mini_kofam/chunk_1/profiles/*.hmm 2>/dev/null | wc -l) HMMs"
echo "  chunk_2: $(ls ${TESTDATA}/mini_kofam/chunk_2/profiles/*.hmm 2>/dev/null | wc -l) HMMs"
