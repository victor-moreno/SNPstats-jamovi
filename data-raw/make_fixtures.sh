#!/usr/bin/env bash
# Generate synthetic genotype fixtures in every format snpImport will read.
#
#   bash data-raw/make_fixtures.sh [tier ...]      default: small
#   bash data-raw/make_fixtures.sh small medium
#   OUT=/somewhere bash data-raw/make_fixtures.sh large
#
# ── Requires PLINK 1.9 ──────────────────────────────────────────────────────
#
#   https://www.cog-genomics.org/plink/1.9/
#
# One tool, on PATH (this script also looks in ~/bin). 1.9 does every format
# conversion here -- .ped/.map, .tped/.tfam, VCF, block-gzipped VCF -- and
# writes the --recode A / --freq / --hardy / --missing oracle files the test
# suite compares against.
#
# plink2 was used for some of this and is no longer needed. It is a different
# tool with different defaults, not a newer version of the same one, and having
# both in the loop meant two installs and two sets of flags for no gain. Three
# places where the 1.9 spelling had to be chosen deliberately to keep the
# fixtures equivalent:
#
#   --recode vcf-iid       plain `vcf` writes FID_IID sample names; plink2's
#                          `--export vcf` wrote the IID alone, and the sample
#                          IDs have to stay the .fam's IIDs
#   --a2-allele            reading a .ped, 1.9 makes A1 the minor allele at
#                          every site; plink2's .bim did not, and a fixture
#                          where it always is tests only half of the allele
#                          orientation handling
#   the generator below    1.9's --dummy gives every variant MAF ~0.5
#
# You do not need plink to run the tests. The `small` tier is committed, and
# that is the only tier tests/testthat uses; this script exists for regenerating
# it and for building the larger tiers, which are gitignored because they run to
# gigabytes. Tests skip themselves when a tier is absent.
#
# ── Provenance ──────────────────────────────────────────────────────────────
#
# Data is completely random, so nothing here is derived from real samples and
# the whole tree is safe to commit or share.
#
# Tiers, and why each format is not generated at every tier:
#
#   small   200 x 500       every format the module reads. Committed; the only
#                           tier tests/testthat uses.
#   medium  2000 x 50000    bed, vcf, tped. The realistic case.
#   large   10000 x 500000  bed only — a .ped at this size is ~20 GB, and
#                           the point of the large tier is to prove that
#                           seeking beats reading, which only needs the .bed.
#
# Only formats the module actually reads are generated. `.gen`/`.sample` was
# produced here while Oxford format was still in scope; it is not, no reader
# exists for it, and no test referenced the files.
#
# Every tier also gets `--recode A` output: the additive (0/1/2) coding of the
# same genotypes, which is the correctness oracle every reader is checked
# against. Cross-format agreement (same SNPs out of .bed and out of .vcf) is a
# second check that needs no oracle at all.
set -euo pipefail

export PATH="$HOME/bin:$PATH"
command -v plink >/dev/null || {
  echo "error: plink (1.9) not on PATH -- https://www.cog-genomics.org/plink/1.9/" >&2
  exit 1; }

HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${OUT:-$HERE/fixtures}"
SEED=20260803

tier_dims() {
  case "$1" in
    small)  echo "200 500" ;;
    medium) echo "2000 50000" ;;
    large)  echo "10000 500000" ;;
    *)      echo "" ;;
  esac
}

# Which text formats to emit per tier. Sizes explode with n*m for the
# sample-major ones, so they are deliberately not generated at every tier.
tier_formats() {
  case "$1" in
    small)  echo "ped tped vcf vcfgz A" ;;
    # vcfgz at medium too: the bgzip test wants matches that fall well past the
    # first ~64 kB block, which small cannot give it. It was missing here, and
    # the test had been running against a medium.vcf.gz made by hand -- so
    # regenerating the tier broke it.
    medium) echo "tped vcf vcfgz A" ;;
    large)  echo "" ;;
  esac
}

make_tier() {
  local tier="$1"
  local dims n m stem fmts
  dims="$(tier_dims "$tier")"
  [ -n "$dims" ] || { echo "unknown tier: $tier" >&2; return 1; }
  n="${dims% *}"; m="${dims#* }"
  stem="$OUT/$tier/$tier"
  fmts="$(tier_formats "$tier")"

  mkdir -p "$OUT/$tier"
  echo "== $tier: $n samples x $m variants -> $OUT/$tier"

  # ── base dataset ──────────────────────────────────────────────────────────
  #
  # The tested tier is written here rather than by --dummy, because --dummy
  # draws every genotype uniformly and so gives *every* variant a minor allele
  # frequency of about 0.5:
  #
  #   plink 1.9 --dummy   MAF  min 0.43  median 0.48  max 0.50
  #   this generator      MAF  min 0.00  median 0.25  max 0.50
  #
  # A panel where nothing is rare cannot test a MAF filter, collapses the
  # difference between an effect-allele frequency and a minor-allele frequency,
  # and gives the HWE test only one kind of site to look at. So the allele
  # frequency spectrum is chosen here instead of inherited from whichever tool
  # generated the data, which is also why this no longer depends on plink2 --
  # its --dummy happened to draw a per-variant frequency, and that accident was
  # load-bearing for six tests.
  #
  # plink still does all the format work: this writes .ped/.map, plink reads it.
  if [ "$tier" = small ]; then
    echo "   generating .ped/.map (chosen MAF spectrum)"
    awk -v m="$m" 'BEGIN{ for (j = 0; j < m; j++) printf "1\tsnp%d\t0\t%d\n", j, j }' \
        > "$stem.map"
    awk -v n="$n" -v m="$m" -v seed="$SEED" -v a2f="$stem.a2" 'BEGIN{
      srand(seed)
      split("A C G T", nt, " ")
      for (j = 0; j < m; j++) {
        a = int(rand()*4)+1; do { b = int(rand()*4)+1 } while (b == a)
        A1[j] = nt[a]; A2[j] = nt[b]
        # which allele plink should put in A2, whatever its frequency
        printf "snp%d\t%s\n", j, nt[b] > a2f
        # Frequency of the FIRST allele, over the whole range rather than
        # capped at 0.5. min(f, 1-f) is then spread over (0.002, 0.5) exactly as
        # before, but the first allele is the minor one only about half the
        # time -- which is what a real .bim looks like, and what stops the
        # effect-allele frequency from being identical to the MAF at every site.
        f[j] = 0.002 + rand()*0.996
      }
      for (i = 0; i < n; i++) {
        # FID flat at 0 on purpose: a cohort with no pedigree structure, which
        # is the shape that exercises the FID-suppression rule end to end. The
        # emitting side is covered by a test that rewrites this .fam with real
        # family IDs.
        line = "0 per" i " 0 0 " (rand()<0.5 ? 1 : 2) " " \
               (rand()<0.05 ? -9 : (rand()<0.5 ? 1 : 2))
        # A few badly genotyped samples, which is what --mind exists to remove.
        # Without them every sample sits within a percent of the same call rate
        # and no threshold can distinguish any of them.
        miss = (rand() < 0.02) ? 0.08 : 0.02
        for (j = 0; j < m; j++) {
          if (rand() < miss) { line = line " 0 0"; continue }
          # Hardy-Weinberg at this variant frequency
          u = rand(); q = f[j]
          if      (u < q*q)             g = A1[j] " " A1[j]
          else if (u < q*q + 2*q*(1-q)) g = A1[j] " " A2[j]
          else                          g = A2[j] " " A2[j]
          line = line " " g
        }
        print line
      }
    }' > "$stem.ped"
    # Reading a .ped, plink assigns A1 = minor allele, and --keep-allele-order
    # does not help: there is no .bim yet whose order it could keep. So A2 is
    # named explicitly, which pins A1 to the other allele whatever its
    # frequency; the generator's second allele is the major one only about half
    # the time.
    #
    # The oracle files below are deliberately written WITHOUT that pinning, so
    # plink reports frequencies for the true minor allele and `.frq`'s MAF
    # column really is a minor allele frequency. The .bim and the .frq therefore
    # name different A1s at about half the sites, which is exactly what the
    # previous plink2-built fixture did.
    #
    # It matters because a .bim where A1 is always the minor allele tests only
    # half of the allele-orientation handling: the effect-allele frequency would
    # equal the MAF at every site, and a weights file naming the major allele
    # would never be exercised. Real files are not like that -- plink2's
    # --make-bed keeps the reference allele in A2 regardless of frequency.
    plink --file "$stem" --a2-allele "$stem.a2" 2 1 \
          --make-bed --out "$stem" >/dev/null 2>&1
    rm -f "$stem.a2"
  else
    # The benchmark tiers are gitignored and exist to measure throughput, where
    # the frequency spectrum does not matter and speed does. --dummy's 0.02 is
    # the per-call missing rate, so the missing code (01 in .bed) is exercised;
    # 1.9 takes the missing *phenotype* rate as a separate argument before the
    # allele-coding modifier, where plink2 took only the genotype rate.
    echo "   plink --dummy"
    plink --dummy "$n" "$m" 0.02 0 acgt \
          --seed "$SEED" \
          --make-bed --out "$stem" >/dev/null 2>&1
    echo "   filling FID, sex and phenotype in .fam"
    awk -v seed="$SEED" 'BEGIN{srand(seed)}
         { $1 = 0;
           $5 = (rand() < 0.5 ? 1 : 2);
           $6 = (rand() < 0.05 ? -9 : (rand() < 0.5 ? 1 : 2));
           print }' "$stem.fam" > "$stem.fam.tmp" && mv "$stem.fam.tmp" "$stem.fam"
  fi

  # ── conversions ───────────────────────────────────────────────────────────
  for f in $fmts; do
    case "$f" in
      ped)
        echo "   -> .ped/.map"
        plink --bfile "$stem" --recode --out "$stem" >/dev/null 2>&1 ;;
      tped)
        # Transposed text is 1.9 only; plink2 dropped it.
        echo "   -> .tped/.tfam"
        plink --bfile "$stem" --recode transpose --out "$stem" >/dev/null 2>&1 ;;
      vcf)
        # vcf-iid, not plain vcf: plain writes FID_IID as the sample name, and
        # the sample IDs have to be the .fam's IIDs so that a VCF import and a
        # .bed import of the same data produce the same IID column.
        echo "   -> .vcf"
        plink --bfile "$stem" --recode vcf-iid --out "$stem" >/dev/null 2>&1 ;;
      vcfgz)
        # The same file block-gzipped, which is what bgzip, htslib and every
        # indexed VCF produce: one independent gzip member per ~64 kB block.
        # DecompressionStream('gzip') decodes only the first, so the reader
        # slices members apart -- and this fixture is what proves it.
        #
        # Written to its own stem and renamed, because --out would otherwise
        # overwrite the plain .vcf produced just above.
        echo "   -> .vcf.gz"
        plink --bfile "$stem" --recode vcf-iid bgz \
              --out "$stem-bgz" >/dev/null 2>&1
        mv "$stem-bgz.vcf.gz" "$stem.vcf.gz"
        rm -f "$stem-bgz.log" "$stem-bgz.nosex" ;;
      A)
        # The oracle: one row per sample, one column per variant, counting the
        # ALT allele. Every reader must reproduce this exactly.
        echo "   -> .raw (--recode A oracle)"
        plink --bfile "$stem" --recode A --out "$stem" >/dev/null 2>&1
        # And plink's own statistics, so the summary table is checked against
        # an outside implementation rather than only against our own oracle:
        #   .frq   MAF        .hwe    genotype counts + exact HWE p
        #   .lmiss per-SNP missingness   .imiss per-sample missingness
        echo "   -> .frq / .hwe / .lmiss / .imiss (summary-table oracles)"
        plink --bfile "$stem" --freq    --out "$stem" >/dev/null 2>&1
        plink --bfile "$stem" --hardy   --out "$stem" >/dev/null 2>&1
        plink --bfile "$stem" --missing --out "$stem" >/dev/null 2>&1 ;;
    esac
  done

  # A short SNP list to select with, and a PGS-Catalog-shaped weights file for
  # the same SNPs — the two ways a selection can be given.
  local k=$(( m < 100 ? m : 100 ))
  echo "   -> selection list and weights file ($k SNPs)"
  awk -v k="$k" 'NR<=k {print $2}' "$stem.bim" > "$OUT/$tier/select.txt"
  { echo "# synthetic PGS weights for testing — not a real score";
    printf "rsID\teffect_allele\tother_allele\teffect_weight\tchr_name\tchr_position\n";
    awk -v k="$k" -v seed="$SEED" 'BEGIN{srand(seed)}
        NR<=k { printf "%s\t%s\t%s\t%.5f\t%s\t%s\n", $2, $5, $6, (rand()*2-1), $1, $4 }' "$stem.bim";
  } > "$OUT/$tier/weights.tsv"

  # A covariate file keyed by IID, for the covariate merge. Deliberately
  # imperfect: it drops a few samples and adds a few that do not exist, so the
  # matched/unmatched reporting has something to report.
  echo "   -> covariates.tsv (with deliberate ID mismatches)"
  { printf "IID\tage\tbmi\tsmoker\n";
    awk -v seed="$SEED" 'BEGIN{srand(seed)}
        NR>3 { if (rand() > 0.03)
                 printf "%s\t%d\t%.1f\t%s\n", $2, 40+int(rand()*40),
                        20+rand()*15, (rand()<0.3 ? "yes" : "no") }' "$stem.fam";
    printf "NOT_A_REAL_SAMPLE_1\t55\t24.0\tno\n";
    printf "NOT_A_REAL_SAMPLE_2\t61\t28.5\tyes\n";
  } > "$OUT/$tier/covariates.tsv"

  echo "   sizes:"
  du -h "$OUT/$tier"/* 2>/dev/null | sort -k2 | sed 's/^/     /'
}

TIERS=("${@:-small}")
for t in "${TIERS[@]}"; do make_tier "$t"; done

echo
echo "done -> $OUT"
