#!/usr/bin/env python3
"""
Set up a full-annotation experiment for the ensembl-genes-nf master pipeline.

Generates a ready-to-run shell script (`run_<assembly>.sh`) containing the
`nextflow run` command with all parameters resolved.

Two modes
---------
Registry mode (recommended on HPC):
    Queries the Ensembl genebuild registry via ensembl-genes Python library to
    resolve assembly metadata, clade settings, stable-ID prefix and stable space.
    Requires the ensembl-genes package to be importable and a settings JSON with
    DB credentials.

    python scripts/setup_experiment.py \\
        --gca             GCA_964261345.1 \\
        --settings-file   ~/genebuild_settings.json \\
        --outdir          /hps/scratch/flicek/ensembl/genebuild/jackt/hetgla

Standalone mode (for testing / off-HPC):
    All metadata supplied manually; no DB connection required.

    python scripts/setup_experiment.py \\
        --gca              GCA_964261345.1 \\
        --assembly-name    mHetGlaV3 \\
        --outdir           /hps/scratch/flicek/ensembl/genebuild/jackt/hetgla \\
        --taxon-id         10181 \\
        --stable-id-prefix ENSHETG \\
        --stable-id-start  1

Registry metadata that cannot be derived automatically (e.g. outdir) must still
be supplied on the command line.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# HPC path constants
# Items marked [CHECK] need verifying on each cluster before the first run.
# ---------------------------------------------------------------------------

# Base directory for IGTR FASTA files (filenames come from clade_settings.json)
IGTR_BASE_DIR = "/hps/nobackup/flicek/ensembl/genebuild/blastdb/ig_tr_proteins"  # [CHECK]

# Selenoprotein FASTA used by the finalise stage
SELENOPROTEIN_FASTA = "/hps/nobackup/flicek/ensembl/genebuild/blastdb/selenoproteins/selenoproteins.fa"  # [CHECK]

# Rfam covariance model file — used by cmsearch in the short_ncrna stage
RFAM_CM = "/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.1/Rfam.cm"  # [CHECK]

# miRNA FASTA — blasted against softmasked genome in short_ncrna
MIRNA_FASTA = "/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/mirBase/all_mirnas.fa"  # [CHECK]

# Projection source genomes and annotations.
# Keys are production_name values from clade_settings.json
# (projection_source_production_name field).
# Each entry needs:
#   fasta — softmasked genome FASTA of the reference species
#   gff3  — canonical gene annotation GFF3 of the reference species
# These files must be pre-generated and available on shared storage.
# Generate the GFF3 with:
#   python ensembl-genes/scripts/export_gff3_from_core.py \
#       --host $GBS5 --port $GBP5 --dbname mus_musculus_core_114_39 \
#       --out mus_musculus.gff3
PROJECTION_SOURCES: dict[str, dict] = {
    "homo_sapiens": {
        "fasta": "/nfs/production/flicek/ensembl/genebuild/projection_sources/homo_sapiens/GRCh38.softmasked.fa",  # [CHECK]
        "gff3":  "/nfs/production/flicek/ensembl/genebuild/projection_sources/homo_sapiens/homo_sapiens.gff3",     # [CHECK]
    },
    "mus_musculus": {
        "fasta": "/nfs/production/flicek/ensembl/genebuild/projection_sources/mus_musculus/GRCm39.softmasked.fa",  # [CHECK]
        "gff3":  "/nfs/production/flicek/ensembl/genebuild/projection_sources/mus_musculus/mus_musculus.gff3",     # [CHECK]
    },
    "danio_rerio": {
        "fasta": "/nfs/production/flicek/ensembl/genebuild/projection_sources/danio_rerio/GRCz11.softmasked.fa",   # [CHECK]
        "gff3":  "/nfs/production/flicek/ensembl/genebuild/projection_sources/danio_rerio/danio_rerio.gff3",       # [CHECK]
    },
    "gallus_gallus": {
        "fasta": "/nfs/production/flicek/ensembl/genebuild/projection_sources/gallus_gallus/bGalGal1.softmasked.fa",  # [CHECK]
        "gff3":  "/nfs/production/flicek/ensembl/genebuild/projection_sources/gallus_gallus/gallus_gallus.gff3",      # [CHECK]
    },
}

# ---------------------------------------------------------------------------
# Augustus species models — closest available trained model per taxon/clade
# ---------------------------------------------------------------------------
AUGUSTUS_SPECIES_MAP: dict[int, str] = {
    # Primates
    9606:  "human", 9544: "human", 9598: "human",
    # Rodents
    10090: "mus_musculus",
    10116: "rattus_norvegicus",
    10181: "mus_musculus",   # Heterocephalus glaber
    9989:  "mus_musculus",   # Rodentia (clade)
    # Other mammals
    9913:  "cow",
    9796:  "horse",
    9823:  "pig",
    9940:  "sheep",
    9615:  "dog",
    40674: "human",          # Mammalia (clade fallback)
    # Birds
    9031:  "chicken", 8782: "chicken",
    # Fish
    7955:  "zebrafish", 7898: "zebrafish",
    # Invertebrates
    7227:  "fly",
    6239:  "caenorhabditis",
    # Plants
    3702:  "arabidopsis", 4530: "rice", 4577: "maize",
    # Fungi
    5141:  "coprinus",
}
AUGUSTUS_DEFAULT = "human"


def get_augustus_species(taxon_id: int | None) -> str:
    if taxon_id and taxon_id in AUGUSTUS_SPECIES_MAP:
        return AUGUSTUS_SPECIES_MAP[taxon_id]
    return AUGUSTUS_DEFAULT


# ---------------------------------------------------------------------------
# Transcriptomic registry lookup
# ---------------------------------------------------------------------------
def query_transcriptomic_registry(
    taxon_id: int,
    settings: dict,
    host: str | None,
    port: str | None,
    min_unique_mapping_pct: float = 60.0,
) -> dict[str, list[str]]:
    """
    Query gb_transcriptomic_registry on the same server as gb_assembly_metadata
    (GBS1 / GBP1) for QC-approved RNA-seq and long-read runs for this taxon.

    Returns:
      {
        "rnaseq":   ["SRR1", "SRR2", ...],   # Illumina, qc_status=ALIGNED
        "longread": ["SRR3", "SRR4", ...],   # ONT / PacBio, qc_status=ALIGNED
        "last_check": "2026-01-08" | None,
      }
    Falls back to empty lists on any error (network, missing table, etc.).

    The min_unique_mapping_pct threshold (default 60 %) further filters runs in
    the align table — runs that mapped poorly to any previous assembly are excluded.
    """
    if not host or not port:
        print(
            "  [INFO] GBS1/GBP1 env vars not set — skipping transcriptome registry query",
            file=sys.stderr,
        )
        return {"rnaseq": [], "longread": [], "last_check": None}

    try:
        import pymysql
        import pymysql.cursors
    except ImportError:
        print(
            "  [INFO] pymysql not installed — skipping transcriptome registry query",
            file=sys.stderr,
        )
        return {"rnaseq": [], "longread": [], "last_check": None}

    try:
        conn = pymysql.connect(
            host=host,
            port=int(port),
            user=settings.get("user_r", settings.get("user", "")),
            password=settings.get("password", ""),
            database="gb_transcriptomic_registry",
            charset="utf8",
            connect_timeout=15,
        )
    except Exception as exc:
        print(f"  [WARN] Could not connect to gb_transcriptomic_registry: {exc}", file=sys.stderr)
        return {"rnaseq": [], "longread": [], "last_check": None}

    rnaseq: list[str]   = []
    longread: list[str] = []
    last_check: str | None = None

    try:
        with conn.cursor(pymysql.cursors.DictCursor) as cur:
            # When was this taxon last checked against ENA?
            cur.execute(
                "SELECT last_check FROM meta WHERE taxon_id = %s LIMIT 1",
                (taxon_id,),
            )
            row = cur.fetchone()
            if row:
                last_check = str(row["last_check"])

            # Fetch ALL qc_status=ALIGNED runs, including their best mapping rate
            # so we can log exactly why any run is excluded.
            cur.execute(
                """
                SELECT DISTINCT
                    r.run_accession,
                    r.platform,
                    MAX(a.uniquely_mapped_reads_percentage) AS best_mapping_pct
                FROM run r
                LEFT JOIN align a ON a.run_id = r.run_id
                WHERE r.taxon_id = %s
                  AND r.qc_status = 'ALIGNED'
                  AND r.read_type = 'RNA-Seq'
                GROUP BY r.run_id, r.run_accession, r.platform
                ORDER BY r.run_id
                """,
                (taxon_id,),
            )
            rejected_runs: list[tuple[str, str, float | None]] = []
            for row in cur.fetchall():
                acc      = row["run_accession"]
                platform = row["platform"]
                best_pct = row["best_mapping_pct"]  # None if no align record yet

                # Accept runs with no align record (not yet mapped, assume OK)
                # Reject runs that mapped poorly to every assembly tried
                if best_pct is not None and best_pct < min_unique_mapping_pct:
                    rejected_runs.append((acc, platform, best_pct))
                    continue

                if platform == "ILLUMINA":
                    rnaseq.append(acc)
                elif platform in ("OXFORD_NANOPORE", "PACBIO_SMRT"):
                    longread.append(acc)

            if rejected_runs:
                print(
                    f"  [REGISTRY] Excluded {len(rejected_runs)} runs "
                    f"(best mapping < {min_unique_mapping_pct}%):",
                    file=sys.stderr,
                )
                for acc, platform, pct in rejected_runs[:10]:
                    print(f"    {acc}  {platform}  {pct:.1f}%", file=sys.stderr)
                if len(rejected_runs) > 10:
                    print(f"    ... and {len(rejected_runs) - 10} more", file=sys.stderr)
    except Exception as exc:
        print(f"  [WARN] Transcriptome registry query failed: {exc}", file=sys.stderr)
    finally:
        conn.close()

    return {"rnaseq": rnaseq, "longread": longread, "last_check": last_check}


# ---------------------------------------------------------------------------
# ENA data lookup (best-effort fallback)
# ---------------------------------------------------------------------------
def _ena_search(taxon_id: int, platform_filter, strategy_filter) -> list[str]:
    """
    Return all unique study_accessions from ENA for the given taxon that match
    the platform and library strategy filters.  Returns [] on failure.
    """
    import urllib.request

    query = f"tax_tree({taxon_id})"   # tax_tree includes subordinate taxa
    fields = "study_accession,instrument_platform,library_strategy,library_source"
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/search"
        f"?display=report&query={query}&domain=read&result=read_run"
        f"&fields={fields}&limit=500&format=json"
    )
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read())
        seen: list[str] = []
        seen_set: set[str] = set()
        for r in data:
            if (
                r.get("instrument_platform") in platform_filter
                and r.get("library_strategy") in strategy_filter
                and r.get("library_source") == "TRANSCRIPTOMIC"
            ):
                bp = r["study_accession"]
                if bp not in seen_set:
                    seen_set.add(bp)
                    seen.append(bp)
        return seen
    except Exception as exc:
        print(f"  [WARN] ENA lookup failed: {exc}", file=sys.stderr)
    return []


def find_rnaseq_bioproject(taxon_id: int) -> str | None:
    """Return comma-separated ENA short-read RNA-seq BioProject(s) or None."""
    hits = _ena_search(
        taxon_id,
        platform_filter={"ILLUMINA"},
        strategy_filter={"RNA-Seq", "EST"},
    )
    return ",".join(hits) if hits else None


def find_longread_bioproject(taxon_id: int) -> str | None:
    """Return comma-separated ENA long-read RNA BioProject(s) or None."""
    hits = _ena_search(
        taxon_id,
        platform_filter={"OXFORD_NANOPORE", "PACBIO_SMRT"},
        strategy_filter={"RNA-Seq"},
    )
    return ",".join(hits) if hits else None


# ---------------------------------------------------------------------------
# Layer priority map
# ---------------------------------------------------------------------------
def build_layer_priorities(stages: dict) -> dict:
    prio: dict[str, int] = {}
    if stages.get("long_read"):   prio["long_read"]          = 0
    if stages.get("targeted"):    prio["best_targeted"]       = 1
    if stages.get("rnaseq"):      prio["rnaseq"]              = 2
    if stages.get("projection"):  prio["projection"]          = 3
    if stages.get("refseq"):      prio["refseq_import"]       = 4
    prio["ab_initio"]             = 5
    if stages.get("genblast"):    prio["genblast_homology"]   = 6
    if stages.get("ncrna"):       prio["short_ncrna"]         = 7
    if stages.get("igtr"):        prio["igtr"]                = 7
    return prio


# ---------------------------------------------------------------------------
# Registry integration (optional — requires ensembl-genes package)
# ---------------------------------------------------------------------------
def _try_import_registry():
    """Return registry helpers or None if ensembl-genes is not installed."""
    try:
        from ensembl.genes.info_from_registry.start_pipeline_from_registry import (
            load_settings,
            get_server_settings_main,
            add_generated_data,
        )
        from ensembl.genes.info_from_registry.assign_species_prefix import (
            get_species_prefix,
        )
        from ensembl.genes.info_from_registry.assign_stable_space import (
            get_stable_space,
        )

        return {
            "load_settings":          load_settings,
            "get_server_settings_main": get_server_settings_main,
            "add_generated_data":     add_generated_data,
            "get_species_prefix":     get_species_prefix,
            "get_stable_space":       get_stable_space,
        }
    except ImportError:
        return None


def resolve_metadata_from_registry(
    gca: str, settings_file: str, outdir: str,
    min_mapping_pct: float = 60.0,
) -> dict:
    """
    Query the Ensembl genebuild registry to resolve all annotation parameters
    for a given GCA accession.

    Returns a flat dict of resolved NF parameters ready for `build_nf_command()`.
    """
    reg = _try_import_registry()
    if reg is None:
        raise RuntimeError(
            "ensembl-genes package not importable. "
            "Run with the full Python environment on the HPC, or use standalone mode "
            "(omit --settings-file and supply --assembly-name / --taxon-id manually)."
        )

    settings = reg["load_settings"](settings_file)

    # Build server_info the same way start_pipeline_from_registry.py main() does:
    # start with the registry DB (GBS1/GBP1 = gb_assembly_metadata), then
    # merge in the main pipeline/core/databases servers.
    # GBS1 / GBP1 are the assembly metadata registry host/port env vars.
    server_info = {
        "registry": {
            "db_host":   os.environ.get("GBS1"),
            "db_user":   settings["user_r"],
            "db_user_w": settings["user"],
            "db_port":   os.environ.get("GBP1"),
            "db_name":   "gb_assembly_metadata",
            "password":  settings["password"],
        }
    }
    # Merge in pipeline_db / core_db / databases keys needed by assign_clade,
    # get_species_prefix, and get_stable_space.
    server_info.update(reg["get_server_settings_main"](settings))

    print(f"Querying registry for {gca}...", file=sys.stderr)
    info = reg["add_generated_data"](server_info, gca, settings)
    # info keys: species_taxon_id, taxon_id, assembly_name, common_name,
    #            assembly_refseq_accession, species_name, assembly_id,
    #            assembly_group, clade, genus_taxon_id, production_name, ...
    #            + clade keys: repbase_library, ig_tr_fasta_file, uniprot_set
    #                          OR protein_file, projection_source_production_name

    taxon_id    = info.get("taxon_id")
    assembly_id = info.get("assembly_id")

    stable_id_prefix = reg["get_species_prefix"](taxon_id, server_info) or ""
    stable_id_start  = reg["get_stable_space"](taxon_id, gca, assembly_id, server_info)
    # stable_id_start is the start integer for this assembly's stable-ID block
    # (get_stable_space returns an int — the stable_space_start value)
    if not isinstance(stable_id_start, int):
        # Function may return the stable_space_id; derive start from it
        stable_id_start = (stable_id_start - 1) * 5_000_000 + 1

    print(f"  Assembly name  : {info['assembly_name']}", file=sys.stderr)
    print(f"  Species        : {info['species_name']}", file=sys.stderr)
    print(f"  Taxon ID       : {taxon_id}", file=sys.stderr)
    print(f"  Clade          : {info.get('clade')}", file=sys.stderr)
    print(f"  Stable-ID pfx  : {stable_id_prefix}", file=sys.stderr)
    print(f"  Stable-ID start: {stable_id_start}", file=sys.stderr)

    # ── Resolve UniProt proteins ──────────────────────────────────────────────
    # Non-vertebrate clades carry an absolute protein_file path.
    # Vertebrate clades carry a uniprot_set name (no path) → use fetch-by-taxon.
    uniprot_fasta    = info.get("protein_file")      # may be None
    uniprot_taxon_id = None
    if not uniprot_fasta:
        # Use the clade-level taxon_id so we get the whole clade's reviewed set
        clade_taxon = info.get("genus_taxon_id") or taxon_id
        uniprot_taxon_id = clade_taxon
        print(
            f"  UniProt source : fetch by taxon {uniprot_taxon_id} "
            f"({info.get('uniprot_set', 'unknown set')})",
            file=sys.stderr,
        )
    else:
        print(f"  UniProt source : {uniprot_fasta}", file=sys.stderr)

    # ── Resolve IGTR proteins (filename → absolute path) ──────────────────────
    igtr_fasta = None
    if info.get("ig_tr_fasta_file"):
        igtr_fasta = str(Path(IGTR_BASE_DIR) / info["ig_tr_fasta_file"])

    # ── Projection source ─────────────────────────────────────────────────────
    projection_source = info.get("projection_source_production_name")
    source_fasta = source_gff3 = None
    if projection_source and projection_source in PROJECTION_SOURCES:
        src = PROJECTION_SOURCES[projection_source]
        source_fasta = src["fasta"]
        source_gff3  = src["gff3"]
        print(f"  Projection     : {projection_source}", file=sys.stderr)
        # Warn if the files don't exist yet — they need to be pre-generated
        for label, path in [("fasta", source_fasta), ("gff3", source_gff3)]:
            if not Path(path).exists():
                print(f"  [WARN] Projection {label} not found: {path}", file=sys.stderr)
                print(f"         Projection stage will be SKIPPED until this file exists.", file=sys.stderr)
    elif projection_source:
        print(
            f"  [WARN] No projection source paths configured for '{projection_source}'. "
            f"Add an entry to PROJECTION_SOURCES in setup_experiment.py.",
            file=sys.stderr,
        )

    # ── Short ncRNA ───────────────────────────────────────────────────────────
    rfam_cm     = RFAM_CM     if Path(RFAM_CM).exists()     else None
    mirna_fasta = MIRNA_FASTA if Path(MIRNA_FASTA).exists() else None
    if not rfam_cm:
        print(f"  [WARN] Rfam CM not found: {RFAM_CM} — short_ncrna stage will be SKIPPED", file=sys.stderr)

    # ── RepeatModeler library ─────────────────────────────────────────────────
    repbase_library = info.get("repbase_library", "vertebrates")

    # ── RefSeq accession ──────────────────────────────────────────────────────
    refseq_accession = info.get("assembly_refseq_accession")

    # ── Augustus species model ────────────────────────────────────────────────
    augustus_species = get_augustus_species(taxon_id)

    # ── Transcriptome registry ────────────────────────────────────────────────
    # Query gb_transcriptomic_registry for QC-approved run accessions.
    # These are preferred over a blind ENA bioproject search because they carry
    # a uniquely_mapped_reads_percentage filter and qc_status=ALIGNED guarantee.
    print("  Querying transcriptome registry...", file=sys.stderr)
    txome = query_transcriptomic_registry(
        taxon_id,
        settings,
        host=os.environ.get("GBS1"),
        port=os.environ.get("GBP1"),
        min_unique_mapping_pct=min_mapping_pct,
    )
    if txome["rnaseq"]:
        print(
            f"  Transcriptome registry: {len(txome['rnaseq'])} QC'd short-read runs "
            f"(last updated: {txome['last_check'] or 'unknown'})",
            file=sys.stderr,
        )
    else:
        print(
            "  Transcriptome registry: no QC'd short-read runs — will fall back to ENA search",
            file=sys.stderr,
        )
    if txome["longread"]:
        print(
            f"  Transcriptome registry: {len(txome['longread'])} QC'd long-read runs",
            file=sys.stderr,
        )

    return {
        "assembly_name":     info["assembly_name"],
        "species_name":      info["species_name"],
        "taxon_id":          taxon_id,
        "clade":             info.get("clade"),
        "stable_id_prefix":  stable_id_prefix,
        "stable_id_start":   stable_id_start,
        "augustus_species":  augustus_species,
        "repbase_library":   repbase_library,
        "refseq_accession":  refseq_accession,
        "uniprot_fasta":     uniprot_fasta,
        "uniprot_taxon_id":  uniprot_taxon_id,
        "igtr_fasta":        igtr_fasta,
        "source_fasta":      source_fasta,
        "source_gff3":       source_gff3,
        "rfam_cm":           rfam_cm,
        "mirna_fasta":       mirna_fasta,
        "outdir":            outdir,
        # Registry-sourced run accessions (comma-sep SRR/ERR IDs), preferred over
        # bioproject-level ENA discovery; None if registry has nothing for this taxon.
        "registry_rnaseq_runs":   ",".join(txome["rnaseq"])   if txome["rnaseq"]   else None,
        "registry_longread_runs": ",".join(txome["longread"]) if txome["longread"] else None,
    }


def resolve_metadata_standalone(args) -> dict:
    """Build the metadata dict from command-line args (no registry access)."""
    if not args.assembly_name:
        sys.exit("--assembly-name is required in standalone mode (no --settings-file).")
    augustus_species = args.augustus_species or get_augustus_species(args.taxon_id)
    rfam_cm     = RFAM_CM     if Path(RFAM_CM).exists()     else None
    mirna_fasta = MIRNA_FASTA if Path(MIRNA_FASTA).exists() else None
    return {
        "assembly_name":          args.assembly_name,
        "species_name":           args.species_name,
        "taxon_id":               args.taxon_id,
        "clade":                  None,
        "stable_id_prefix":       args.stable_id_prefix or "",
        "stable_id_start":        args.stable_id_start,
        "augustus_species":       augustus_species,
        "repbase_library":        args.repbase_library,
        "refseq_accession":       args.refseq_accession,
        "uniprot_fasta":          args.uniprot_fasta,
        "uniprot_taxon_id":       args.uniprot_taxon_id,
        "igtr_fasta":             None,
        "source_fasta":           args.source_fasta,
        "source_gff3":            args.source_gff3,
        "rfam_cm":                rfam_cm,
        "mirna_fasta":            mirna_fasta,
        "outdir":                 args.outdir.rstrip("/"),
        # No registry access in standalone mode
        "registry_rnaseq_runs":   None,
        "registry_longread_runs": None,
    }


# ---------------------------------------------------------------------------
# Build the nextflow run command
# ---------------------------------------------------------------------------
def build_nf_command(
    *,
    gca: str,
    meta: dict,
    repo_dir: str,
    nf_work_root: str,
    profile: str,
    rnaseq_nf_arg: str,
    longread_nf_arg: str,
    has_rnaseq: bool,
    has_longread: bool,
    has_genblast: bool,
    skip_repeatmodeler: bool,
) -> tuple[str, dict]:
    """
    Assemble the `nextflow run` command string and the Prefect parameter dict.

    Returns (nf_cmd_string, prefect_params_dict).
    """
    has_projection = bool(meta.get("source_fasta") and meta.get("source_gff3"))
    has_ncrna      = bool(meta.get("rfam_cm"))

    stages = {
        "long_read":  has_longread,
        "rnaseq":     has_rnaseq,
        "projection": has_projection,
        "genblast":   has_genblast,
        "refseq":     bool(meta.get("refseq_accession")),
        "ncrna":      has_ncrna,
        "igtr":       bool(meta.get("igtr_fasta")),
    }
    layer_prios_json = json.dumps(build_layer_priorities(stages))

    nf_lines = [
        f"nextflow run {repo_dir}",
        f"  --assembly_accession        {gca}",
        f"  --assembly_name             {meta['assembly_name']}",
        f"  --outdir                    {meta['outdir']}",
        f"  --nf_work_root              {nf_work_root}",
        f"  --ab_initio_species         {meta['augustus_species']}",
        f"  --repbase_library           {meta['repbase_library']}",
        f"  --layer_priorities          '{layer_prios_json}'",
        f"  --stable_id_prefix          '{meta['stable_id_prefix']}'",
        f"  --stable_id_start           {meta['stable_id_start']}",
    ]

    # RNA-seq
    if rnaseq_nf_arg:
        nf_lines.append(f"  {rnaseq_nf_arg}")

    # Long-read
    if longread_nf_arg:
        nf_lines.append(f"  {longread_nf_arg}")

    # Projection
    if has_projection:
        nf_lines.append(f"  --source_fasta              {meta['source_fasta']}")
        nf_lines.append(f"  --source_gff3               {meta['source_gff3']}")

    # Short ncRNA
    if meta.get("rfam_cm"):
        nf_lines.append(f"  --rfam_cm                   {meta['rfam_cm']}")
    if meta.get("mirna_fasta"):
        nf_lines.append(f"  --mirna_fasta               {meta['mirna_fasta']}")

    # UniProt / protein homology
    if meta.get("uniprot_fasta"):
        nf_lines.append(f"  --uniprot_fasta             {meta['uniprot_fasta']}")
    elif meta.get("uniprot_taxon_id"):
        nf_lines.append(f"  --uniprot_taxon_id          {meta['uniprot_taxon_id']}")

    # IGTR
    if meta.get("igtr_fasta"):
        nf_lines.append(f"  --igtr_proteins             {meta['igtr_fasta']}")

    # RefSeq
    if meta.get("refseq_accession"):
        nf_lines.append(f"  --assembly_refseq_accession {meta['refseq_accession']}")

    # Selenoproteins
    nf_lines.append(f"  --selenoprotein_fasta       {SELENOPROTEIN_FASTA}")

    # RepeatModeler
    if skip_repeatmodeler:
        nf_lines.append("  --skip_repeatmodeler        true")

    # Species name (for core DB loading)
    if meta.get("species_name"):
        nf_lines.append(f"  --species_name              '{meta['species_name']}'")

    nf_lines.append(f"  -profile                    {profile}")
    nf_lines.append("  -resume")

    nf_cmd = " \\\n".join(nf_lines)

    prefect_params = {
        "_migration_note": (
            "This JSON will drive a Prefect flow once the eHive → Prefect migration "
            "is complete. Each key maps to a Prefect task input."
        ),
        "assembly_accession":        gca,
        "assembly_name":             meta["assembly_name"],
        "outdir":                    meta["outdir"],
        "taxon_id":                  meta.get("taxon_id"),
        "clade":                     meta.get("clade"),
        "ab_initio_species":         meta["augustus_species"],
        "repbase_library":           meta["repbase_library"],
        "layer_priorities":          build_layer_priorities(stages),
        "rnaseq_source":             rnaseq_nf_arg or None,
        "longread_source":           longread_nf_arg or None,
        "source_fasta":              meta.get("source_fasta"),
        "source_gff3":               meta.get("source_gff3"),
        "rfam_cm":                   meta.get("rfam_cm"),
        "mirna_fasta":               meta.get("mirna_fasta"),
        "uniprot_fasta":             meta.get("uniprot_fasta"),
        "uniprot_taxon_id":          meta.get("uniprot_taxon_id"),
        "igtr_proteins":             meta.get("igtr_fasta"),
        "assembly_refseq_accession": meta.get("refseq_accession"),
        "stable_id_prefix":          meta["stable_id_prefix"],
        "stable_id_start":           meta["stable_id_start"],
    }

    return nf_cmd, prefect_params


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ── Required ──────────────────────────────────────────────────────────────
    parser.add_argument("--gca",     required=True, help="GCA accession (e.g. GCA_964261345.1)")
    parser.add_argument("--outdir",  required=True, help="Root output directory on HPC scratch")

    # ── Registry mode ─────────────────────────────────────────────────────────
    registry_grp = parser.add_argument_group(
        "Registry mode",
        "Supply --settings-file to query the Ensembl genebuild registry. "
        "All other metadata will be resolved automatically.",
    )
    registry_grp.add_argument(
        "--settings-file", default=None,
        help="JSON settings file with DB credentials (same format as ensembl-genes pipeline scripts)",
    )
    registry_grp.add_argument(
        "--min-mapping-pct", type=float, default=60.0,
        help="Minimum STAR uniquely_mapped_reads_percentage for a run to pass from the "
             "transcriptome registry (default: 60.0)",
    )

    # ── Standalone overrides ───────────────────────────────────────────────────
    standalone_grp = parser.add_argument_group(
        "Standalone / override",
        "Used in standalone mode (no --settings-file) or to override registry values.",
    )
    standalone_grp.add_argument("--assembly-name",    default=None, help="Assembly name (no spaces)")
    standalone_grp.add_argument("--taxon-id",         type=int, default=None, help="NCBI taxon ID")
    standalone_grp.add_argument("--species-name",     default=None, help="Scientific name")
    standalone_grp.add_argument("--stable-id-prefix", default=None, help="Stable-ID prefix (e.g. ENSHETG)")
    standalone_grp.add_argument("--stable-id-start",  type=int, default=1, help="Stable-ID range start")
    standalone_grp.add_argument("--repbase-library",  default="vertebrates", help="RepeatMasker library")
    standalone_grp.add_argument("--refseq-accession", default=None, help="GCF accession for RefSeq import")
    standalone_grp.add_argument("--augustus-species", default=None, help="Override Augustus species model")
    standalone_grp.add_argument("--skip-repeatmodeler", action="store_true")
    # Projection source (standalone or registry override)
    standalone_grp.add_argument("--source-fasta", default=None,
                                help="Projection source genome FASTA (overrides registry lookup)")
    standalone_grp.add_argument("--source-gff3",  default=None,
                                help="Projection source annotation GFF3 (overrides registry lookup)")

    # ── RNA-seq (mutually exclusive) ──────────────────────────────────────────
    g_rna = parser.add_mutually_exclusive_group()
    g_rna.add_argument("--sample-sheet",         default=None, help="Local CSV sample sheet")
    g_rna.add_argument("--rnaseq-bioproject",    default=None, help="ENA BioProject accession(s), comma-separated")
    g_rna.add_argument("--rnaseq-run-accessions",default=None, help="Comma-separated SRR/ERR accessions")

    # ── Long-read RNA (mutually exclusive) ────────────────────────────────────
    g_lr = parser.add_mutually_exclusive_group()
    g_lr.add_argument("--long-read-bioproject", default=None,
                      help="ENA BioProject accession(s) for long-read RNA (comma-separated)")
    g_lr.add_argument("--long-read-csv",        default=None,
                      help="Pre-built long-read CSV sample sheet (tab-separated, ENA format)")

    # ── UniProt (mutually exclusive, standalone override) ─────────────────────
    g_prot = parser.add_mutually_exclusive_group()
    g_prot.add_argument("--uniprot-fasta",    default=None, help="Local UniProt FASTA (overrides registry)")
    g_prot.add_argument("--uniprot-taxon-id", type=int, default=None,
                        help="Taxon ID → auto-fetch from UniProt (overrides registry)")

    # ── Infrastructure ────────────────────────────────────────────────────────
    parser.add_argument("--nf-work-root", default=None,  help="NF work dir (default: outdir/.nf_work)")
    parser.add_argument("--repo-dir",     default=".",   help="Path to ensembl-genes-nf repo")
    parser.add_argument("--profile",      default="cluster", help="Nextflow profile")
    parser.add_argument("--output-script",default=None,  help="Output .sh path (default: run_<assembly>.sh)")

    args = parser.parse_args()

    outdir       = args.outdir.rstrip("/")
    nf_work_root = args.nf_work_root or f"{outdir}/.nf_work"
    repo_dir     = str(Path(args.repo_dir).resolve())

    # ── Resolve metadata ───────────────────────────────────────────────────────
    if args.settings_file:
        meta = resolve_metadata_from_registry(args.gca, args.settings_file, outdir,
                                               min_mapping_pct=args.min_mapping_pct)
        # Allow CLI overrides on top of registry values
        if args.assembly_name:   meta["assembly_name"]   = args.assembly_name
        if args.taxon_id:        meta["taxon_id"]        = args.taxon_id
        if args.species_name:    meta["species_name"]    = args.species_name
        if args.stable_id_prefix:meta["stable_id_prefix"]= args.stable_id_prefix
        if args.stable_id_start != 1: meta["stable_id_start"] = args.stable_id_start
        if args.repbase_library != "vertebrates": meta["repbase_library"] = args.repbase_library
        if args.refseq_accession:meta["refseq_accession"]= args.refseq_accession
        if args.augustus_species:meta["augustus_species"] = args.augustus_species
        if args.uniprot_fasta:
            meta["uniprot_fasta"]    = args.uniprot_fasta
            meta["uniprot_taxon_id"] = None
        if args.uniprot_taxon_id:
            meta["uniprot_taxon_id"] = args.uniprot_taxon_id
            meta["uniprot_fasta"]    = None
        # Projection source overrides
        if args.source_fasta:
            meta["source_fasta"] = args.source_fasta
        if args.source_gff3:
            meta["source_gff3"] = args.source_gff3
    else:
        meta = resolve_metadata_standalone(args)

    skip_repeatmodeler = args.skip_repeatmodeler

    assembly_name = meta["assembly_name"]
    script_path   = args.output_script or f"run_{assembly_name}.sh"

    print(f"Augustus species model : {meta['augustus_species']}", file=sys.stderr)

    # ── Resolve RNA-seq ────────────────────────────────────────────────────────
    # Priority: explicit CLI arg > transcriptome registry > ENA bioproject search
    rnaseq_nf_arg = ""
    has_rnaseq    = False

    if args.sample_sheet:
        rnaseq_nf_arg = f"--sample_sheet              {args.sample_sheet}"
        has_rnaseq    = True
    elif args.rnaseq_bioproject:
        rnaseq_nf_arg = f"--rnaseq_bioproject         {args.rnaseq_bioproject}"
        has_rnaseq    = True
    elif args.rnaseq_run_accessions:
        rnaseq_nf_arg = f"--rnaseq_run_accessions     {args.rnaseq_run_accessions}"
        has_rnaseq    = True
    elif meta.get("registry_rnaseq_runs"):
        # Use QC-approved run list from genebuild transcriptome registry
        runs = meta["registry_rnaseq_runs"]
        n    = len(runs.split(","))
        print(f"  Using {n} QC'd short-read runs from genebuild registry", file=sys.stderr)
        rnaseq_nf_arg = f"--rnaseq_run_accessions     {runs}"
        has_rnaseq    = True
    elif meta.get("taxon_id"):
        print(
            f"  No registry RNA-seq — querying ENA for taxon {meta['taxon_id']}...",
            file=sys.stderr,
        )
        bioproject = find_rnaseq_bioproject(meta["taxon_id"])
        if bioproject:
            n_projects = len(bioproject.split(","))
            print(f"  Found {n_projects} ENA project(s): {bioproject}", file=sys.stderr)
            rnaseq_nf_arg = f"--rnaseq_bioproject         {bioproject}"
            has_rnaseq    = True
        else:
            print("  [WARN] No short-read RNA-seq found — skipping RNA-seq stage.", file=sys.stderr)

    # ── Resolve long-read RNA ──────────────────────────────────────────────────
    # Priority: explicit CLI arg > transcriptome registry > ENA bioproject search
    longread_nf_arg = ""
    has_longread    = False

    if args.long_read_csv:
        longread_nf_arg = f"--long_read_csv             {args.long_read_csv}"
        has_longread    = True
    elif args.long_read_bioproject:
        longread_nf_arg = f"--long_read_bioproject      {args.long_read_bioproject}"
        has_longread    = True
    elif meta.get("registry_longread_runs"):
        runs = meta["registry_longread_runs"]
        n    = len(runs.split(","))
        print(f"  Using {n} QC'd long-read runs from genebuild registry", file=sys.stderr)
        # Pass individual run accessions so fetch_reads_from_ena.py downloads
        # each run separately with the correct platform filter
        longread_nf_arg = f"--long_read_bioproject      {runs}"
        has_longread    = True
    elif meta.get("taxon_id"):
        print(
            f"  No registry long-read — querying ENA for taxon {meta['taxon_id']}...",
            file=sys.stderr,
        )
        lr_bioproject = find_longread_bioproject(meta["taxon_id"])
        if lr_bioproject:
            n_projects = len(lr_bioproject.split(","))
            print(f"  Found {n_projects} ENA long-read project(s): {lr_bioproject}", file=sys.stderr)
            longread_nf_arg = f"--long_read_bioproject      {lr_bioproject}"
            has_longread    = True
        else:
            print("  [INFO] No long-read RNA found — skipping long-read stage.", file=sys.stderr)

    has_genblast = bool(meta.get("uniprot_fasta") or meta.get("uniprot_taxon_id"))

    # ── Build command ──────────────────────────────────────────────────────────
    nf_cmd, prefect_params = build_nf_command(
        gca                = args.gca,
        meta               = meta,
        repo_dir           = repo_dir,
        nf_work_root       = nf_work_root,
        profile            = args.profile,
        rnaseq_nf_arg      = rnaseq_nf_arg,
        longread_nf_arg    = longread_nf_arg,
        has_rnaseq         = has_rnaseq,
        has_longread       = has_longread,
        has_genblast       = has_genblast,
        skip_repeatmodeler = skip_repeatmodeler,
    )

    # ── Write run script ───────────────────────────────────────────────────────
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    script_content = f"""#!/bin/bash
# ensembl-genes-nf full annotation
# Generated : {now}
# Assembly  : {args.gca} ({assembly_name})
# Taxon     : {meta.get('taxon_id') or 'not specified'}
# Clade     : {meta.get('clade') or 'not specified'}
#
# Launch inside tmux on the cluster:
#   tmux new -s {assembly_name.lower()}
#   salloc --ntasks=1 --mem=4G --time=7-00:00:00 --partition=standard
#   bash {script_path}

set -euo pipefail

echo "Starting annotation for {args.gca} at $(date)"
mkdir -p {outdir}

{nf_cmd}

EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    echo "Pipeline completed successfully at $(date)"
    echo "Final annotation: {outdir}/finalise_geneset/finalise_geneset/*.canonical.gff3"
else
    echo "Pipeline FAILED with exit code $EXIT_CODE at $(date)" >&2
    exit $EXIT_CODE
fi
"""

    Path(script_path).write_text(script_content)
    os.chmod(script_path, 0o755)

    # ── Write Prefect migration stub ───────────────────────────────────────────
    prefect_path = f"prefect_params_{assembly_name}.json"
    Path(prefect_path).write_text(json.dumps(prefect_params, indent=2))

    # ── Summary ────────────────────────────────────────────────────────────────
    has_projection = bool(meta.get("source_fasta") and meta.get("source_gff3"))
    has_ncrna      = bool(meta.get("rfam_cm"))
    stage_flags = {
        "long_read":  has_longread,
        "rnaseq":     has_rnaseq,
        "projection": has_projection,
        "genblast":   has_genblast,
        "refseq":     bool(meta.get("refseq_accession")),
        "ncrna":      has_ncrna,
        "igtr":       bool(meta.get("igtr_fasta")),
    }
    layer_prios_json = json.dumps(build_layer_priorities(stage_flags))

    def _tick(flag: bool) -> str:
        return "✓" if flag else "✗"

    print(f"\n{'='*65}")
    print(f"  Experiment : {args.gca} ({assembly_name})")
    if meta.get("clade"):
        print(f"  Clade      : {meta['clade']}")
    print(f"{'='*65}")
    print(f"  Run script     : {script_path}")
    print(f"  Prefect params : {prefect_path}")
    print(f"  Outdir         : {outdir}")
    print(f"\n  Active stages:")
    print(f"    {_tick(True)}  Repeat masking   : RepeatModeler + RepeatMasker ({meta['repbase_library']})")
    print(f"    {_tick(True)}  Ab initio        : Augustus '{meta['augustus_species']}'")
    if has_rnaseq:
        print(f"    {_tick(True)}  RNA-seq          : {rnaseq_nf_arg.strip()}")
    else:
        print(f"    {_tick(False)}  RNA-seq          : not found / not provided")
    if has_longread:
        print(f"    {_tick(True)}  Long-read RNA    : {longread_nf_arg.strip()}")
    else:
        print(f"    {_tick(False)}  Long-read RNA    : not found / not provided")
    if has_projection:
        print(f"    {_tick(True)}  Projection       : {meta.get('source_fasta', '')}")
    else:
        print(f"    {_tick(False)}  Projection       : source files not available")
    if has_genblast:
        src = meta.get("uniprot_fasta") or f"fetch taxon {meta.get('uniprot_taxon_id')}"
        print(f"    {_tick(True)}  Protein homology : {src}")
    else:
        print(f"    {_tick(False)}  Protein homology : no UniProt source")
    if meta.get("igtr_fasta"):
        print(f"    {_tick(True)}  IG/TR            : {meta['igtr_fasta']}")
    else:
        print(f"    {_tick(False)}  IG/TR            : not configured")
    if has_ncrna:
        print(f"    {_tick(True)}  Short ncRNA      : Rfam + miRNA")
    else:
        print(f"    {_tick(False)}  Short ncRNA      : Rfam CM not found")
    if meta.get("refseq_accession"):
        print(f"    {_tick(True)}  RefSeq import    : {meta['refseq_accession']}")
    else:
        print(f"    {_tick(False)}  RefSeq import    : no RefSeq accession")
    print(f"\n  Layer priorities : {layer_prios_json}")
    print(f"\n  Stable ID prefix : {meta['stable_id_prefix']}")
    print(f"  Stable ID start  : {meta['stable_id_start']}")
    print(f"\n  To run:")
    print(f"    tmux new -s {assembly_name.lower()}")
    print(f"    salloc --ntasks=1 --mem=4G --time=7-00:00:00 --partition=standard")
    print(f"    bash {script_path}")
    print(f"{'='*65}\n")


if __name__ == "__main__":
    main()
