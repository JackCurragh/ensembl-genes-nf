# Annotation Pipeline

Translon annotation workflow with two complementary approaches:

## Overview

This pipeline provides two strategies for translon annotation, each serving different purposes:

### 1. Global Annotation (Mainstream)
**Purpose**: Leverage published tools for standard genome-wide annotation

**Approach**: Sample clustering → published tools → consensus
- Clusters all samples globally into N meta-samples
- Runs established tools (ORFquant, RiboCode, Ribotish, etc.)
- Each tool operates on full genome context
- Generates consensus annotation per cluster
- **Scale**: Run each tool ~10 times (N clusters)

**Best for**: Standard translon discovery, tool comparison, benchmark datasets

### 2. Per-Gene Annotation (Custom)
**Purpose**: Find gene-specific, condition-dependent translation patterns

**Approach**: Gene-level clustering → custom methods → profiles
- Clusters samples separately for each gene based on translation profiles
- Generates multiple annotations per gene (one per cluster)
- Uses custom profile-based and statistical methods
- Detects condition-specific translons
- **Scale**: Process ~20,000 genes, each with multiple clusters

**Best for**: Differential translation analysis, condition-specific discoveries, heterogeneous datasets

## Quick Start

### Run Both Branches (Default)

```bash
nextflow run pipelines/annotation/main.nf \
  --mode both \
  --gtf reference.gtf \
  --fasta reference.fa \
  --riboseq_results_dir results/riboseq/ \
  --count_matrix results/riboseq/count_matrix.parquet \
  --merged_bam results/riboseq/merged.bam \
  --sample_offsets results/riboseq/offsets/
```

### Global Annotation Only

```bash
nextflow run pipelines/annotation/main.nf \
  --mode global \
  --gtf reference.gtf \
  --fasta reference.fa \
  --riboseq_results_dir results/riboseq/ \
  --translon_callers orfquant,ribocode,ribotish
```

### Per-Gene Annotation Only

```bash
nextflow run pipelines/annotation/main.nf \
  --mode pergene \
  --gtf reference.gtf \
  --fasta reference.fa \
  --count_matrix results/riboseq/count_matrix.parquet \
  --merged_bam results/riboseq/merged.bam \
  --sample_offsets results/riboseq/offsets/
```

## Parameters

### Common Parameters
- `--mode`: Pipeline mode (`global`, `pergene`, or `both`)
- `--gtf`: Reference GTF annotation file (required)
- `--fasta`: Reference genome FASTA file (required)
- `--rnaseq_results_dir`: Optional RNA-Seq results (enhances both branches)

### Global Annotation Parameters
- `--riboseq_results_dir`: Ribo-Seq pipeline results directory
- `--similarity_method`: Method for clustering (`correlation`, `cosine`, `euclidean`)
- `--similarity_threshold`: Global clustering threshold (default: 0.8)
- `--translon_callers`: Comma-separated list (default: `orfquant,ribocode,ribotish`)
- `--min_caller_agreement`: Minimum callers for consensus (default: 2)

### Per-Gene Annotation Parameters
- `--count_matrix`: Unique RPFs × samples matrix (parquet format)
- `--merged_bam`: All samples merged BAM file
- `--sample_offsets`: Per-sample offset files
- `--pergene_similarity_method`: Similarity method for per-gene clustering
- `--pergene_similarity_threshold`: Per-gene clustering threshold (default: 0.8)
- `--min_cluster_size`: Minimum samples per gene-cluster (default: 3)
- `--min_periodicity_score`: Minimum 3-nt periodicity (default: 0.5)
- `--min_coverage`: Minimum reads per translon (default: 10)
- `--p_value_threshold`: Statistical significance (default: 0.05)

## Outputs

```
results/annotation/
├── global/
│   ├── cluster_assignments.tsv         # Sample → cluster mapping
│   ├── merged_bams/                    # Merged BAMs per cluster
│   ├── caller_predictions/             # Individual tool outputs
│   │   ├── orfquant/
│   │   ├── ribocode/
│   │   └── ribotish/
│   ├── consensus/
│   │   ├── cluster_1.gtf               # Final annotation for cluster 1
│   │   ├── cluster_2.gtf               # Final annotation for cluster 2
│   │   └── annotation_stats.json
│   └── reports/
│       └── global_annotation_report.html
│
└── pergene/
    ├── gene_cluster_assignments.tsv    # Gene → cluster → samples mapping
    ├── cluster_profiles/               # Translation profiles per gene-cluster
    ├── candidate_translons/            # Raw predictions per gene
    ├── scored_translons/               # Filtered predictions with scores
    ├── final_annotations/
    │   ├── gene_cluster_annotations.tsv
    │   └── pergene_annotation.gtf      # All annotations combined
    └── reports/
        └── pergene_annotation_report.html
```

## Workflow Details

### Global Annotation Workflow

1. **Load Ribo-Seq Data**: Import BAMs, offsets, metadata
2. **Compute Profiles**: Generate translation profiles across all samples
3. **Cluster Samples**: Global hierarchical clustering by similarity
4. **Merge Data**: Combine BAMs and offsets within each cluster
5. **Filter Transcripts**: Optional RNA-Seq-based filtering
6. **Call Translons**: Run published tools on each cluster
7. **Standardize**: Convert all outputs to common format
8. **Consensus**: Generate final annotation per cluster

### Per-Gene Annotation Workflow

1. **Load Matrix**: Import unique RPF count matrix
2. **Cluster Per-Gene**: Hierarchical clustering for each gene independently
3. **Extract BAMs**: Slice gene-specific reads for each cluster
4. **Generate Profiles**: Create cluster-specific translation profiles
5. **Profile Detection**: Custom periodicity, enrichment, frame analysis
6. **Statistical Tests**: Compare clusters for differential translation
7. **Score Translons**: Combine evidence, filter by confidence
8. **Annotate**: Generate per-gene, per-cluster annotations

## Implementation Status

⚠️ **Both branches are currently templated at the subworkflow level.**

- Subworkflow structure defined
- Module placeholders created
- Ready for implementation

See individual subworkflow files for detailed TODO steps.

## Design Rationale

**Why two branches?**

- **Global**: Published tools aren't designed for per-gene operation. They need full genome context for offset calculation, parameter learning, etc.
- **Per-Gene**: Finds biology that global methods miss - condition-specific translation that varies gene-by-gene

**When to use which?**

- Use **Global** for: Standard annotation, method benchmarking, tool comparison
- Use **Per-Gene** for: Differential translation, heterogeneous samples, condition discovery
- Use **Both** for: Comprehensive analysis - global gives standard annotation, per-gene finds special cases
