# Test fixtures

Synthetic genotype data in every format `snpImport` will read, plus the files
that drive SNP selection and the covariate merge.

**Everything here is random.** Nothing is derived from real samples, so the
whole tree is safe to commit and share.

## Generating

```bash
bash data-raw/make_fixtures.sh                 # small (default)
bash data-raw/make_fixtures.sh small medium
bash data-raw/make_fixtures.sh large
```

Needs [`plink` 1.9](https://www.cog-genomics.org/plink/1.9/) on `PATH`; the
script adds `~/bin`. Output is deterministic — same seed, same files — so a
regenerated tier is byte-identical, which is checked by regenerating into a
scratch directory and comparing.

The `small` tier is written by the script itself as `.ped`/`.map` and handed to
plink, rather than coming from `--dummy`: `--dummy` draws every genotype
uniformly, so every variant ends up with a minor allele frequency near 0.5, and
a panel where nothing is rare cannot test a MAF filter. The benchmark tiers do
use `--dummy`, where throughput is the point and the frequency spectrum is not.

| tier | samples × variants | formats | size |
|---|---|---|---|
| small | 200 × 500 | all of them | ~1.5 MB, **committed** |
| medium | 2 000 × 50 000 | bed, tped, vcf | ~25 MB, gitignored |
| large | 10 000 × 500 000 | bed only | ~1.2 GB, gitignored |

The text formats are not generated at every tier on purpose: a `.ped` at the
large tier would be roughly 20 GB, and the large tier exists to prove that
seeking beats reading, which needs only the `.bed`.

## What each tier contains

```
<tier>.bed/.bim/.fam      the base dataset
<tier>.ped/.map           sample-major text
<tier>.tped/.tfam         variant-major text
<tier>.vcf                VCF
<tier>.vcf.gz             bgzipped VCF (BGZF, several members)
<tier>.raw                --recode A: the correctness oracle
select.txt                100 variant IDs, one per line
weights.tsv               the same 100 SNPs as a PGS-Catalog weights file
covariates.tsv            IID + age/bmi/smoker, with deliberate ID mismatches
```

`--dummy` is called with a 2 % missing rate and `acgt`, so the fixtures
exercise the `.bed` missing code (`01`) and real allele letters rather than
`1`/`2`. The `.fam` is post-processed to carry a realistic sex code and a
1/2 case-control phenotype with ~5 % `-9` missing, because `snpImport` reads
both out of it.

`covariates.tsv` deliberately drops ~3 % of the genotyped samples and adds two
IDs that do not exist. The covariate merge has to report both
directions, and a fixture where everything matches would not test that.

## What checks these against plink

The oracle files are not decoration — the test suite compares against them cell
by cell, which is what makes the readers falsifiable:

| test file | oracle | what it proves |
|---|---|---|
| `test-oracle.R` | `.raw` (`--recode A`) | every decoded genotype, cell by cell |
| `test-golden-plink.R` | `.frq`, `.hwe`, `.lmiss`, `.imiss` | MAF, HWE and missingness match plink's own |
| `test-formats.R` | none needed | `.bed`, `.ped`, `.tped` and VCF agree with each other |

Run them with `bash tests/run_tests.sh`, or `bash tests/run-in-docker.sh` where
the R library is not available.

One thing the `.raw` comparison taught us, worth keeping in mind for the
importer: **`--recode A` counts the minor allele, which is A1 for only about
half the variants** in random data. Comparing against it needs a per-variant
orientation, not a uniform `2 - x` flip. That is the same allele-orientation
trap the weights-file path has to handle, and it is why `test-oracle.R` reads
the counted allele out of the column header rather than assuming it.
