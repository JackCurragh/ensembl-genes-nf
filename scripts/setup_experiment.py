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
# HPC path constants — flag [CHECK] items before running on a new cluster
# ---------------------------------------------------------------------------

# Base directory for IGTR FASTA files referenced by filename in clade_settings.json
IGTR_BASE_DIR = (
    "/hps/nobackup/flicek/ensembl/genebuild/blastdb/ig_tr_proteins"  # [CHECK]
)

# RepeatModeler library base (only used if --skip-repeatmodeler is NOT set
# and no pre-built library is found; RepeatModeler will build from scratch)
REPEATMODELER_LIB_BASE = (
    "/hps/nobackup/flicek/ensembl/genebuild/repeatmodeler_libraries"  # [CHECK]
)

# Selenoprotein FASTA used by the finalise stage
SELENOPROTEIN_FASTA = (
    "/hps/nobackup/flicek/ensembl/genebuild/blastdb/selenoproteins/selenoproteins.fa"  # [CHECK]
)

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
# ENA RNA-seq lookup (best-effort — only used in standalone mode)
# ---------------------------------------------------------------------------
def find_rnaseq_bioproject(taxon_id: int) -> str | None:
    import urllib.request

    query  = f"tax_eq({taxon_id})"
    fields = "study_accession,instrument_platform,library_strategy"
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/search"
        f"?display=report&query={query}&domain=read&result=read_run"
        f"&fields={fields}&limit=10&format=json"
    )
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read())
        illumina = [
            r for r in data
            if r.get("instrument_platform") == "ILLUMINA"
            and r.get("library_strategy") in ("RNA-Seq", "EST")
        ]
        if illumina:
            return illumina[0]["study_accession"]
    except Exception as exc:
        print(f"  [WARN] ENA lookup failed: {exc}", file=sys.stderr)
    return None


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
    gca: str, settings_file: str, outdir: str
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

    # ── RepeatModeler library ─────────────────────────────────────────────────
    repbase_library = info.get("repbase_library", "vertebrates")

    # ── RefSeq accession ──────────────────────────────────────────────────────
    refseq_accession = info.get("assembly_refseq_accession")

    # ── Augustus species model ────────────────────────────────────────────────
    augustus_species = get_augustus_species(taxon_id)

    return {
        "assembly_name":       info["assembly_name"],
        "species_name":        info["species_name"],
        "taxon_id":            taxon_id,
        "clade":               info.get("clade"),
        "stable_id_prefix":    stable_id_prefix,
        "stable_id_start":     stable_id_start,
        "augustus_species":    augustus_species,
        "repbase_library":     repbase_library,
        "refseq_accession":    refseq_accession,
        "uniprot_fasta":       uniprot_fasta,
        "uniprot_taxon_id":    uniprot_taxon_id,
        "igtr_fasta":          igtr_fasta,
        "projection_source":   info.get("projection_source_production_name"),
        "outdir":              outdir,
    }


def resolve_metadata_standalone(args) -> dict:
    """Build the metadata dict from command-line args (no registry access)."""
    if not args.assembly_name:
        sys.exit("--assembly-name is required in standalone mode (no --settings-file).")
    augustus_species = args.augustus_species or get_augustus_species(args.taxon_id)
    return {
        "assembly_name":    args.assembly_name,
        "species_name":     args.species_name,
        "taxon_id":         args.taxon_id,
        "clade":            None,
        "stable_id_prefix": args.stable_id_prefix or "",
        "stable_id_start":  args.stable_id_start,
        "augustus_species": augustus_species,
        "repbase_library":  args.repbase_library,
        "refseq_accession": args.refseq_accession,
        "uniprot_fasta":    args.uniprot_fasta,
        "uniprot_taxon_id": args.uniprot_taxon_id,
        "igtr_fasta":       None,
        "projection_source": None,
        "outdir":           args.outdir.rstrip("/"),
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
    has_rnaseq: bool,
    has_genblast: bool,
    skip_repeatmodeler: bool,
) -> tuple[str, dict]:
    """
    Assemble the `nextflow run` command string and the Prefect parameter dict.

    Returns (nf_cmd_string, prefect_params_dict).
    """
    stages = {
        "rnaseq":   has_rnaseq,
        "genblast": has_genblast,
        "refseq":   bool(meta["refseq_accession"]),
        "igtr":     bool(meta.get("igtr_fasta")),
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

    # UniProt / protein homology
    if meta.get("uniprot_fasta"):
        nf_lines.append(f"  --uniprot_fasta             {meta['uniprot_fasta']}")
    elif meta.get("uniprot_taxon_id"):
        nf_lines.append(f"  --uniprot_taxon_id          {meta['uniprot_taxon_id']}")

    # IGTR
    if meta.get("igtr_fasta"):
        nf_lines.append(f"  --igtr_proteins             {meta['igtr_fasta']}")

    # RefSeq
    if meta["refseq_accession"]:
        nf_lines.append(f"  --assembly_refseq_accession {meta['refseq_accession']}")

    # RepeatModeler
    if skip_repeatmodeler:
        nf_lines.append("  --skip_repeatmodeler        true")

    # Species name (for core DB loading)
    if meta.get("species_name"):
        nf_lines.append(f"  --species_name              '{meta['species_name']}'")

    # Selenoproteins
    if Path(SELENOPROTEIN_FASTA).exists() or True:  # always include; NF checks existence
        nf_lines.append(f"  --selenoprotein_fasta       {SELENOPROTEIN_FASTA}")

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
        "uniprot_fasta":             meta.get("uniprot_fasta"),
        "uniprot_taxon_id":          meta.get("uniprot_taxon_id"),
        "igtr_proteins":             meta.get("igtr_fasta"),
        "assembly_refseq_accession": meta["refseq_accession"],
        "stable_id_prefix":          meta["stable_id_prefix"],
        "stable_id_start":           meta["stable_id_start"],
        "projection_source":         meta.get("projection_source"),
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

    # ── RNA-seq (mutually exclusive) ──────────────────────────────────────────
    g_rna = parser.add_mutually_exclusive_group()
    g_rna.add_argument("--sample-sheet",         default=None, help="Local CSV sample sheet")
    g_rna.add_argument("--rnaseq-bioproject",    default=None, help="ENA BioProject accession")
    g_rna.add_argument("--rnaseq-run-accessions",default=None, help="Comma-separated SRR/ERR accessions")

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
        meta = resolve_metadata_from_registry(args.gca, args.settings_file, outdir)
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
    else:
        meta = resolve_metadata_standalone(args)

    skip_repeatmodeler = args.skip_repeatmodeler

    assembly_name = meta["assembly_name"]
    script_path   = args.output_script or f"run_{assembly_name}.sh"

    print(f"Augustus species model : {meta['augustus_species']}", file=sys.stderr)

    # ── Resolve RNA-seq ────────────────────────────────────────────────────────
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
    elif meta.get("taxon_id"):
        print(
            f"No RNA-seq specified — querying ENA for taxon {meta['taxon_id']}...",
            file=sys.stderr,
        )
        bioproject = find_rnaseq_bioproject(meta["taxon_id"])
        if bioproject:
            print(f"  Found: {bioproject}", file=sys.stderr)
            rnaseq_nf_arg = f"--rnaseq_bioproject         {bioproject}"
            has_rnaseq    = True
        else:
            print("  [WARN] No RNA-seq found in ENA — running ab initio only.", file=sys.stderr)

    has_genblast = bool(meta.get("uniprot_fasta") or meta.get("uniprot_taxon_id"))

    # ── Build command ──────────────────────────────────────────────────────────
    nf_cmd, prefect_params = build_nf_command(
        gca                = args.gca,
        meta               = meta,
        repo_dir           = repo_dir,
        nf_work_root       = nf_work_root,
        profile            = args.profile,
        rnaseq_nf_arg      = rnaseq_nf_arg,
        has_rnaseq         = has_rnaseq,
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
    layer_prios_json = json.dumps(build_layer_priorities({
        "rnaseq":   has_rnaseq,
        "genblast": has_genblast,
        "refseq":   bool(meta["refseq_accession"]),
        "igtr":     bool(meta.get("igtr_fasta")),
    }))

    print(f"\n{'='*65}")
    print(f"  Experiment : {args.gca} ({assembly_name})")
    if meta.get("clade"):
        print(f"  Clade      : {meta['clade']}")
    print(f"{'='*65}")
    print(f"  Run script     : {script_path}")
    print(f"  Prefect params : {prefect_path}")
    print(f"  Outdir         : {outdir}")
    print(f"\n  Active stages:")
    print(f"    Ab initio     : Augustus '{meta['augustus_species']}' (always)")
    if has_rnaseq:
        print(f"    RNA-seq       : {rnaseq_nf_arg.strip()}")
    if has_genblast:
        src = meta.get("uniprot_fasta") or f"fetch taxon {meta.get('uniprot_taxon_id')}"
        print(f"    Protein homol : {src}")
    if meta.get("igtr_fasta"):
        print(f"    IGTR          : {meta['igtr_fasta']}")
    if meta["refseq_accession"]:
        print(f"    RefSeq import : {meta['refseq_accession']}")
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
