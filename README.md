# genome-annotation

A Nextflow DSL2 pipeline for functional annotation of prokaryotic and microbial genomes. Starting from genome FASTA files, it predicts protein-coding sequences and annotates them against Pfam, KEGG (KOfam), and CAZy (dbCAN) databases. Optionally reconstructs genome-scale metabolic models using GapSeq and/or CarVeME.

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D)](https://www.nextflow.io/)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Pipeline overview

```
Genome FASTA(s)
      │
      ▼
  SEQSTATS          ← count bases/sequences; branch small (<100 kb) vs large
      │
      ▼
  PYRODIGAL         ← gene prediction (meta mode for small genomes, single for large)
      │
      ├──▶ HMMER_HMMSEARCH   ── Pfam domain annotation     ──▶ <sample>_<db>.domtbl.gz
      │
      ├──▶ RUNDBCAN           ── CAZyme annotation          ──▶ overview / HMM / substrate / DIAMOND TSVs
      │
      ├──▶ KOFAMSCAN          ── KEGG ortholog assignment   ──▶ <sample>_kofam.txt
      │
      └──▶ CARVEME_CARVE*     ── metabolic model (SBML)     ──▶ <sample>.xml

Genome FASTA ──▶ GAPSEQ_DOALL* ── metabolic model (RDS + SBML)
```

\* Optional, enabled by `--run_carveme` / `--run_gapseq`.

## Quick start

### 1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
```

### 2. Prepare a samplesheet

`samplesheet.csv`:
```csv
sample,fastx
genome_A,/path/to/genome_A.fa
genome_B,/path/to/genome_B.fa.gz
```

Each row is one genome. The `fastx` column accepts uncompressed or gzip-compressed FASTA.

### 3. Run

Database paths can be passed on the command line or via a params file (recommended — see [Configuration](#configuration)):

```bash
nextflow run main.nf \
    -profile docker \
    --samplesheet samplesheet.csv \
    --outdir results \
    --databases.pfam.hmm /path/to/Pfam-A.hmm \
    --databases.dbcan.files.db /path/to/dbcan_db/ \
    --databases.kofam_chunked.chunk_1.files.ko_list /path/to/ko_list \
    --databases.kofam_chunked.chunk_1.files.profiles /path/to/profiles/
```

## Databases

| Database | Where to download | Param |
|----------|-------------------|-------|
| Pfam-A HMM | [InterPro FTP](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/) — `Pfam-A.hmm.gz` | `params.databases.pfam.hmm` |
| dbCAN v14+ | [dbCAN2 downloads](https://bcb.unl.edu/dbCAN2/download/Databases/) — full directory | `params.databases.dbcan.files.db` |
| KOfam | [KEGG FTP](https://www.genome.jp/ftp/db/kofam/) — `profiles.tar.gz` + `ko_list` | `params.databases.kofam_chunked.*` |

KOfam can be split into multiple chunks to parallelise annotation. Define as many `chunk_N` blocks as needed under `params.databases.kofam_chunked` (see [Configuration](#configuration)).

### Metabolic modelling databases (optional)

| Tool | Database | Param |
|------|----------|-------|
| GapSeq | [GapSeq data](https://github.com/jotech/gapseq) — run `gapseq fetch dat` | `params.databases.gapseq.path` |
| CarVeME | BIGG universe model (bundled in container) | — |

## Configuration

The recommended approach is a params file:

```nextflow
// my_params.config
params {
    samplesheet = '/path/to/samplesheet.csv'
    outdir      = '/path/to/results'

    databases {
        pfam {
            hmm        = '/data/Pfam-A.hmm'
            n_profiles = 24736   // update if using a non-standard Pfam release
        }
        dbcan {
            files { db = '/data/dbcan/' }
        }
        kofam_chunked {
            chunk_1 {
                files {
                    ko_list  = '/data/kofam/ko_list'
                    profiles = '/data/kofam/profiles/'
                }
            }
            // add chunk_2, chunk_3 … to parallelise across database partitions
        }
    }

    run_gapseq  = false
    run_carveme = false
}
```

Run with:
```bash
nextflow run main.nf -profile docker -params-file my_params.config
```

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--samplesheet` | — | Path to input CSV samplesheet |
| `--outdir` | `results` | Output directory |
| `--prodigal_meta` | `false` | Force Pyrodigal meta mode for all genomes (useful for MAGs) |
| `--pfam_chunksize` | `2 MB` | Protein FASTA chunk size for Pfam HMMER jobs |
| `--kofam_chunksize` | `2 MB` | Protein FASTA chunk size for KOfam jobs |
| `--dbcan_chunksize` | `1 MB` | Protein FASTA chunk size for dbCAN jobs |
| `--run_gapseq` | `false` | Run GapSeq metabolic model reconstruction |
| `--run_carveme` | `false` | Run CarVeME metabolic model reconstruction |
| `--gapseq_medium` | — | Growth medium CSV for GapSeq (auto-predicted if unset) |
| `--carveme_mediadb` | — | CarVeME media database file |

## Execution profiles

| Profile | Description |
|---------|-------------|
| `docker` | Docker (AMD64; applies `--platform linux/amd64` for Apple Silicon compatibility) |
| `singularity` | Singularity/Apptainer |
| `test` | Minimal test run using bundled mini databases (run `tests/scripts/make_mini_databases.sh` once first) |
| `local` | Resource limits for a workstation |
| `codon` | SLURM profile for the EMBL-EBI CODON HPC cluster |

```bash
# Test run
nextflow run main.nf -profile test,docker --outdir test_results
```

## Output structure

```
results/
└── <sample_id>/
    ├── pyrodigal/
    │   ├── <sample>.faa.gz               # predicted proteins
    │   ├── <sample>.fna.gz               # predicted CDS nucleotides
    │   └── <sample>.gff.gz               # gene coordinates (GFF3)
    ├── pfam/
    │   └── concatenated/
    │       └── <sample>_<db>.domtbl.gz   # HMMER domain table
    ├── dbcan/
    │   └── concatenated/
    │       ├── <sample>_dbcan_overview.tsv
    │       ├── <sample>_dbcan_cazyme.tsv
    │       ├── <sample>_dbcan_substrate.tsv
    │       └── <sample>_dbcan_diamond.tsv
    ├── kofam/
    │   └── concatenated/
    │       └── <sample>_kofam.txt
    ├── gapseq/                           # present if --run_gapseq
    │   ├── <sample>.RDS                  # gap-filled model (R object)
    │   ├── <sample>.xml                  # gap-filled SBML model
    │   ├── <sample>-draft.RDS
    │   ├── <sample>-draft.xml
    │   └── <sample>-medium.csv           # predicted growth medium
    └── carveme/                          # present if --run_carveme
        └── <sample>.xml                  # CarVeME SBML model
```

## Testing

```bash
# One-time: build mini test databases from full databases in ~/Downloads/
bash tests/scripts/make_mini_databases.sh

# Run all tests (all module + subworkflow + pipeline tests)
nf-test test --profile +docker

# Run only the end-to-end pipeline test
nf-test test tests/default.nf.test --profile +docker

# Run a specific module test
nf-test test modules/local/seqstats/tests/main.nf.test --profile +docker
```

## Utility scripts

### `scripts/inject_kofam_ga_thresholds.py`

Injects KOfam score thresholds from a `ko_list` file into HMM profiles as HMMER3 GA (Gathering Threshold) lines. This enables use of `--cut_ga` with `hmmsearch` for automatic hit filtering, without a separate post-processing step.

Both GA fields (full-sequence and domain) are set to the KOfam threshold value. No external dependencies — requires only Python 3.

```bash
python scripts/inject_kofam_ga_thresholds.py \
    --ko-list /path/to/ko_list \
    --profiles-dir /path/to/profiles/ \
    --output-dir /path/to/profiles_with_ga/
```

Profiles for KOs with no defined threshold (marked `-` in `ko_list`) are copied unchanged. Use `--in-place` to modify the profiles directory directly rather than writing to a new location.

```
options:
  --ko-list FILE          KOfam ko_list TSV file
  --profiles-dir DIR      Directory of individual .hmm profile files
  --output-dir DIR        Write modified profiles here (mutually exclusive with --in-place)
  --in-place              Overwrite profiles in-place
  --missing-threshold     {skip,warn,error}  behaviour for profiles with no ko_list entry (default: skip)
```

Run the unit tests with:
```bash
pytest scripts/tests/test_inject_kofam_ga_thresholds.py -v
```

## Credits

Written by Tim Rozday ([@timrozday-mgnify](https://github.com/timrozday-mgnify)).

Built on the [nf-core](https://nf-co.re) framework and community modules. Please cite the nf-core publication if you use this pipeline:

> Ewels PA, Peltzer A, Fillinger S, et al. **The nf-core framework for community-curated bioinformatics pipelines.** *Nat Biotechnol.* 2020;38:276–278. doi:[10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)

Tool citations: [`CITATIONS.md`](CITATIONS.md).
