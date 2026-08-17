# QOV — virus and prophage identification in human gut metagenome assemblies

QOV takes a metagenome assembly and returns the viral and prophage sequences it contains, each with
an admission score and, optionally, a predicted bacterial host.

Everything the pipeline needs is inside one Apptainer image: no conda, no environment activation,
no external database downloads. The image is ~5.5 GB and carries CheckV, geNomad, MetaPhaPred,
DeepVirFinder, PHROGs, the viral-vs-cellular contrast reference and the admission model, all baked
in. The single exception is the optional Kraken2 host database (19 GB), which the image detects
automatically if you bind it and skips if you do not.

## Download the prebuilt image

The image is distributed prebuilt, so nothing needs compiling:

```bash
curl -L -o QOVirome.sif https://datacloud.helsinki.fi/index.php/s/JfFByXAwgYbBkQ9/download
```

(Browser link: <https://datacloud.helsinki.fi/index.php/s/JfFByXAwgYbBkQ9>)

The optional Kraken2 host database is a separate download and is only needed for host prediction.

## Usage

```bash
apptainer run --cleanenv QOVirome.sif --in contigs.fna --outdir results --threads 8
```

That is the whole invocation. No binds, no database staging, no environment variables.

With predicted bacterial host (bind the optional Kraken2 database — host prediction switches itself
on when it is present):

```bash
apptainer run --cleanenv \
  --bind /path/to/kraken2/HRGMv2:/opt/gutbusters/db/kraken2/HRGMv2:ro \
  QOVirome.sif --in contigs.fna --outdir results --threads 8
```

| flag | default | meaning |
|---|---|---|
| `--in` | *(required)* | assembly FASTA |
| `--outdir` | `results` | output directory |
| `--threads` | `4` | worker threads |
| `--min-len` | `2000` | minimum contig length in bp |
| `--score` | `0.9` | phage-score parameter |
| `--prefix` | `busted` | output file prefix |
| `--keeptmp` | off | retain intermediate files |

Force host prediction on or off explicitly with `--env GB_KRAKEN=1` / `--env GB_KRAKEN=0`.
Override the admission threshold with `--env GB_STRICTNESS=sensitive|balanced|strict`, or set an
exact cut with `--env GB_VSCORE_MIN=0.5140`.

### Outputs

- `results/<prefix>_viruses.fna` — admitted viral and prophage sequences. Excised prophage regions
  keep their parent contig name with a `_ckv_<region>` suffix.
- `results/<prefix>_table.tsv` — one row per admitted sequence:

  `ID, length, dvf, mpp, h_rank, h_name, h_score, viral_genes, host_genes, viral_frac, bact_frac, viral_score`

  `viral_score` is the admission score in [0, 1]; the `h_*` columns are populated only when the
  Kraken2 database is bound.
- `results/logs/` — one log per pipeline stage.

A virus-free sample is a valid result: both files are written empty and the run exits 0.

## Admission thresholds

| level | cut | precision | recall |
|---|---|---|---|
| sensitive | 0.4040 | 79.71% | 77.16% |
| balanced *(default)* | 0.5140 | 87.14% | 66.96% |
| strict | 0.6560 | 94.61% | 54.25% |

Measured out-of-sample on real gut metagenome contigs. Precision and recall on reference-genome
fragments are considerably higher; see the manuscript for why that number should not be quoted on
its own.

## How it works

1. **Gather.** Permissive tools claim candidate viral sequence, and CheckV excises prophage regions
   from bacterial contigs. No single claimer is adequate on its own — the best reaches 44% recall —
   but their union reaches 99%.
2. **Describe.** Each candidate is described by 176 features that combine *measured evidence*
   (CheckV gene census, PHROGs profile hits, viral-vs-cellular protein bitscore contrast, geNomad
   marker scores) with *learned prediction* (DeepVirFinder, MetaPhaPred, tetranucleotide
   composition).
3. **Decide.** One supervised model (extremely randomised trees, 800 estimators) scores every
   candidate, and the score is thresholded at the chosen strictness.

Scoring at runtime walks the frozen model in pure Python, without scikit-learn, so a library
upgrade cannot silently move a score.

## Also available as a hosted service

The same pipeline runs at <https://questomics.app> with no installation at all.

## Repository layout

```
Apptainer.selfcontained.def   image definition — the authoritative record of every tool and database version
pipelines/                    the three executable pipeline stages, in order
scripts/                      feature extraction and the model scorer
entrypoint/gutbusters         stage dispatcher
run.sh                        production wrapper (S3 staging, database checks); not needed for local use
benchmarks/                   the derived tables behind every figure in the manuscript
```

## Building the image

Requires the databases listed in the `%files` section of the definition. Build on a host with
~26 GB of free scratch space:

```bash
singularity build QOVirome.sif Apptainer.selfcontained.def
```

## Citation

See `CITATION.cff`. The manuscript is in preparation; this repository will be archived with a
Zenodo DOI on submission.


## The distributed model

The admission model is an ExtraTrees classifier serialised to JSON and baked into the image, so
scoring is deterministic and needs no Python model object at runtime.

| | |
|---|---|
| features | 158 |
| estimator | `ExtraTreesClassifier(n_estimators=800, max_features=0.2, criterion='entropy', random_state=20260811)` |
| training set | 14,916 contigs from 45 donors, labelled by homology to the Unified Human Gut Virome |
| mandatory signals | `gn_virus`, `gn_chr`, `gn_plasmid`, `mpp` — a contig missing any of these is dropped rather than imputed |
| imputable signals | the remaining 154, each with a frozen training median |

It is a **single fit**, not an ensemble: every number reported in the manuscript is reproducible
from this one artifact.

## Licence

See `LICENSE`. Third-party tools bundled in the image retain their own licences.
