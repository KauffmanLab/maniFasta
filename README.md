# maniFasta

**maniFasta builds reproducible protein reference databases for metaproteomics from heterogeneous source types.**

maniFasta accepts heterogeneous inputs - including species lists, protein accession lists, and user-supplied FASTA files - and combines their proteins into a harmonized reference database.

Each build produces two core outputs: a protein **FASTA** with standardized, source-aware headers and a per-protein provenance **mani**fest. maniFasta also packages the configuration, inputs, scripts, and logs used for the build and generates a draft methods summary for review.

[AllOralsDB](https://www.homd.org/ftp/AllOralsDB/) is the motivating application of maniFasta. It was developed as a taxonomically comprehensive reference database for salivary and oral metaproteomics and will continue to be empirically refined, with new versions posted as they are developed.

## Get started

### Running maniFasta on Google Colab

The interactive guided maniFasta notebook supports source selection, user uploads, validation, taxonomy-output controls, live build output, and download of the complete run. No prepackaged sources are selected by default. The example button selects one source per module to demonstrate the workflow. Colab storage is temporary, so be sure to download the completed ZIP of your build before closing the runtime.

**Latest stable version** (`v2026.220.1007`): [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/KauffmanLab/maniFasta/blob/v2026.220.1007/maniFasta_on_GoogleColab.ipynb), or the **current version** (`main`): [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/KauffmanLab/maniFasta/blob/main/maniFasta_on_GoogleColab.ipynb)

### Running maniFasta at command line

maniFasta requires Linux, Bash, Python 3.10 or newer, and standard Unix command-line tools. The [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) is required for NCBI genome and proteome downloads.

```bash
git clone https://github.com/KauffmanLab/maniFasta.git
cd maniFasta

cp 00.setup/000.maniFasta.local.env.template \
   00.setup/000.maniFasta.local.env
chmod 600 00.setup/000.maniFasta.local.env
```

Then:

1. Set `mainDIR`, `run_label`, and `db_prefix` in `00.setup/000.maniFasta.config`.
2. To use prepackaged sources, set `enabled=TRUE` for the desired rows in `00.setup/000.maniFasta.source_registry.tsv`. All sources are disabled by default.
3. To use your own accession lists, protein lists, or FASTA files, place them in `00.setup/`, add a source-registry row that references each file, and set that row to `enabled=TRUE`. Field guidance and examples are provided in the registry.
4. Add the path to `datasets` and optional NCBI credentials to `00.setup/000.maniFasta.local.env` when required by the selected sources.

Run directly:

```bash
cd 01.scripts
export Config_file="$(pwd)/../00.setup/000.maniFasta.config"
bash runBuild_db.sh
```

Or submit with Slurm after reviewing the scheduler and resource settings in the config:

```bash
cd 01.scripts
bash 01.submit_runBuild_db.sh
```

## Choose sources

The source registry defines the database composition. One row represents one source, and enabled rows may use any of four modules. Each module reflects a different input type that maniFasta can handle:

| Module | Input |
|---|---|
| `mod_I` | List(s) of proteomes defined by organism name, taxid, genome or proteome accession, or UniProt proteome ID |
| `mod_II` | List(s) of individual protein accessions from NCBI, UniProt, UniParc, or PDB |
| `mod_III` | Protein FASTA file(s), with optional accompanying per-protein metadata |
| `mod_IV` | Supported collections: human UniProt, HOMD PROKKA proteomes, or cRAP contaminants |

Prepackaged inputs are located in `99.prepackaged_inputs/`. Detailed fields and options for prepackaged and custom sources are documented directly in the config and source registry.

## Outputs

Each run creates `YYYYMMDD_HHMMSS_<run_label>/` containing:

- `<run_label>.fasta`: the harmonized protein reference database
- `<run_label>.manifest.tsv`: per-protein source and taxonomy provenance
- `build_summary.tsv`: sequence totals and build status
- `<run_label>.database_methods.md`: a draft methods summary for review
- optional lineage-enriched metadata and SVG/HTML taxonomy sunbursts
- `run_data/`, `run_setup/`, `run_scripts/`, and `run_logs/`: the data and documentation used for the build
- `maniFasta_software_provenance.tsv`: the exact maniFasta Git tag and commit used

## AllOralsDB

AllOralsDB applies maniFasta to salivary and oral metaproteomics. It includes human proteins, oral bacteria and archaea, fungi and other microeukaryotes, viruses and virus-like elements (including Obelisks), selected dairy proteins, and common contaminants.

Four <I>e</I>HOMD flavors balance taxonomic breadth against database size:

| Flavor | HOMD selection |
|---|---|
| A | All oral-associated assemblies |
| S | One representative assembly per oral species |
| G | One representative assembly per oral genus |
| F | One representative assembly per oral family |

Current AllOralsDB releases and their associated documentation are available from [<I>e</I>HOMD](https://www.homd.org/ftp/AllOralsDB/).

## Important behavior

- maniFasta retains all successfully retrieved or staged records unless a source-specific filter excludes them; it does not deduplicate the final FASTA.
- Custom FASTA residue alphabets and terminal stop symbols are not automatically sanitized.
- Missing accessions may be skipped and reported or treated as fatal through `allow_missing_accessions` in the config.
- PDB retrieval currently supports entry IDs and polymer-entity IDs, not chain-form IDs; the affected AllergenOnline v24 rows are documented in the registry.
- Optional lineage and taxonomy reporting may require substantial memory for multi-million-protein databases.

## Citation and contributors

When using maniFasta or AllOralsDB, please cite our [preprint](https://www.biorxiv.org/content/10.64898/2026.08.08.739415v1), as well as the original references for any prepackaged sources included in your database.

maniFasta was developed by [Christopher Handelmann](https://scholar.google.com/citations?user=1WPxKpsAAAAJ&hl=en&oi=ao), [Ashley K. Miles](https://akmiles-code.github.io/), and [Kathryn M. Kauffman](https://kauffmanlab.org/)
