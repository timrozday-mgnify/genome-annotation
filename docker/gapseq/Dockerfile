FROM mambaorg/micromamba:2.0.5 AS builder

USER root

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install all gapseq runtime dependencies, then strip build-only artifacts
# in the same layer so they never appear in the runtime image's COPY layer.
RUN micromamba install -y -n base --copy \
        -c conda-forge \
        -c bioconda \
        bash \
        r-base \
        gawk \
        sed \
        grep \
        bc \
        git \
        coreutils \
        wget \
        openssl \
        glpk \
        hmmer \
        blast \
        "diamond>=2.1.9" \
        mmseqs2 \
        pyrodigal \
        libsbml \
        r-data.table \
        r-stringr \
        r-stringi \
        r-getopt \
        r-r.utils \
        "r-cobrar>=0.2.2" \
        r-biocmanager \
        bioconductor-biostrings \
        r-jsonlite \
        r-renv \
        r-rcurl \
        r-httr \
    && micromamba clean --all --yes \
    && rm -rf /opt/conda/pkgs \
    # GCC cross-compile toolchain (~917 MB) — used only during conda package
    # compilation, not needed to run R/BLAST/Diamond/MMseqs2/HMMER at runtime
    && rm -rf \
        /opt/conda/x86_64-conda-linux-gnu \
        /opt/conda/libexec/gcc \
        /opt/conda/lib/gcc \
    && rm -f /opt/conda/bin/x86_64-conda-linux-gnu-* \
    # Static archives — already linked into binaries at build time
    # maxdepth 1 for top-level .a, recurse into ncbi-blast+ for its static archives
    && find /opt/conda/lib -maxdepth 1 -name '*.a' -delete \
    && find /opt/conda/lib -maxdepth 1 -name '*.a.*' -delete \
    && find /opt/conda/lib/ncbi-blast+ -name '*.a' -delete \
    # Sanitizer DSOs — only used with -fsanitize= instrumented builds
    && rm -f \
        /opt/conda/lib/libasan.so* \
        /opt/conda/lib/libtsan.so* \
        /opt/conda/lib/liblsan.so* \
        /opt/conda/lib/libhwasan.so* \
        /opt/conda/lib/libubsan.so* \
    # entrez-direct tools — not used by gapseq
    && rm -f \
        /opt/conda/bin/efetch \
        /opt/conda/bin/esearch \
        /opt/conda/bin/elink \
        /opt/conda/bin/epost \
        /opt/conda/bin/efilt \
        /opt/conda/bin/transmute.Linux \
        /opt/conda/bin/rchive.Linux \
        /opt/conda/bin/xtract.Linux \
        /opt/conda/bin/datatool \
        /opt/conda/bin/nghttpx \
        /opt/conda/bin/project_tree_builder \
    # Header files — only needed to compile C/C++ code
    && rm -rf /opt/conda/include \
    # share: remove locale strings, documentation, man/info pages, easel C sources
    && rm -rf \
        /opt/conda/share/locale \
        /opt/conda/share/doc \
        /opt/conda/share/man \
        /opt/conda/share/info \
        /opt/conda/share/gir-1.0 \
        /opt/conda/share/easel \
        /opt/conda/share/git-gui \
        /opt/conda/share/glib-2.0/codegen \
    # Perl runtime — only used by legacy BLAST Perl scripts, not called by gapseq
    && rm -rf /opt/conda/lib/perl5 /opt/conda/bin/perl /opt/conda/bin/perl5.32.0 \
    && rm -f /opt/conda/bin/legacy_blast.pl /opt/conda/bin/update_blastdb.pl \
       /opt/conda/bin/aggregate_profile.pl /opt/conda/bin/profile2mat.pl \
    # Python bytecode caches and test suite — not needed to run pyrodigal
    && find /opt/conda -name '*.pyc' -delete \
    && find /opt/conda -name '__pycache__' -type d -empty -delete \
    && rm -rf /opt/conda/lib/python3.13/test \
    # Miscellaneous non-runtime data
    && rm -rf \
        /opt/conda/man \
        /opt/conda/fonts \
        /opt/conda/conda-meta

# ── runtime stage ────────────────────────────────────────────────────────────
FROM mambaorg/micromamba:2.0.5

LABEL org.opencontainers.image.title="gapseq"
LABEL org.opencontainers.image.version="2.0.1"
LABEL org.opencontainers.image.description="Informed prediction of bacterial metabolic pathways and reconstruction of accurate metabolic models"
LABEL org.opencontainers.image.source="https://github.com/jotech/gapseq"
LABEL org.opencontainers.image.licenses="GPL-3.0"

USER root

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Copy the already-cleaned conda environment from the builder stage
COPY --from=builder /opt/conda /opt/conda

# procps (provides `ps`) is used by some gapseq scripts
RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

# Clone gapseq v2.0.1 from source, then:
#   - Remove files not needed at runtime (~115 MB)
#   - Strip bundled sequence DB (Bacteria/Archaea user seqs must be provided via -D)
#   - Patch prepare_batch_alignments.R so user/ sequences use the external -D path
RUN git clone --depth 1 --branch v2.0.1 https://github.com/jotech/gapseq /opt/gapseq \
    && chmod +x /opt/gapseq/gapseq \
    && rm -rf \
        /opt/gapseq/.git \
        /opt/gapseq/toy \
        /opt/gapseq/docs \
        /opt/gapseq/unit \
        /opt/gapseq/dat/seq/Bacteria \
        /opt/gapseq/dat/seq/Archaea \
    && mkdir -p /opt/gapseq/dat/seq/Bacteria/user /opt/gapseq/dat/seq/Archaea/user \
    && sed -i \
        's|fileDT\[, sdir := ifelse(src=="user",paste0(script.dir,"/../dat/seq"),seqdb)\]|fileDT[, sdir := seqdb]|' \
        /opt/gapseq/src/prepare_batch_alignments.R

ENV PATH="/opt/conda/bin:/opt/gapseq:${PATH}"
