# GWAS Summary Stats File Format Definition

| COL # | COL NAME | TYPE | DESCRIPTION | GENESIS | RVTESTS |
| :---: | --- | --- | --- | --- | --- |
| 0 | `ID` | STR | Variant ID | `variant.id`<sup>1</sup> |  |
| 1 | `CHR` | INT | Chromosome of variant | `chr` | `CHROM` |
| 2 | `POS` | INT | Position of variant | `pos` | `POS` |
| 3 | `REF` | STR | Reference allele | `REF(0)` from info | `REF` |
| 4 | `ALT` | STR | Alternate allele | `ALT(1)` from info | `ALT` |
| 5 | `ALT_AF` | FLOAT | Frequency of the alternate allele | `freq`<sup>1</sup> | `ALT_Frq` |
| 6 | `MAF` | FLOAT | Frequency of the minor allele | `MAF` from info | `MAF` from info |
| 7 | `POP_MAF` | FLOAT | Population frequency of the minor allele | `AFR`, `AMR`, `EAS`, `EUR` or `SAS` from 1000G legend | `AFR`, `AMR`, `EAS`, `EUR` or `SAS` from 1000G legend |
| 8 | `SOURCE` | STR | Source of variant genotypes (`obs` or `imp`) | `Genotyped` from info | `Genotyped` from info |
| 9 | `IMP_QUAL` | FLOAT | Imputation R<sup>2</sup> | `Rsq` from info | `Rsq` from info |
| 10 | `N` | INT | Number of individuals contributing to analysis | `n.obs` | `N_INFORMATIVE` |
| 11 | `ALT_EFFECT` | FLOAT | Effect size of alternate allele | `Est`<sup>1</sup> | `ALT_EFFSIZE` |
| 12 | `SE` | FLOAT | Standard error | `Est.SE` | `SQRT_V_STAT`<sup>1</sup> |
| 13 | `P-VALUE` | FLOAT | P-value of association | `Score.pval` | `PVALUE` |

<sup>1</sup> Some modification necessary

For discussion:
- Additional fields?
- ID format - `[CHR]:[POSITION]:[REF]:[ALT]` vs. RS ID vs. Oxford vs. our own identifier
- Do we need both `ALT_AF` and `MAF`? Or do we just include `MAF` with a standard that the `ALT` allele is always the minor allele?
- Standard for specifying indel alleles (e.g., `GA|G` vs `A|-`) and position
  - I think we'll need to do the latter, because it's possible to unambiguously convert from the 1st to the 2nd, but not the other way (unless the extra allele in the 1st is invariant - not sure if this is true or not). On the other hand, making this conversion would result in a loss of information.
  - Position depends on allele convention we use (I think)


Unused GENESIS fields:
- MAC: the minor allele count
- Score: the value of the score function
- Score.SE: the estimated standard error of the score
- Score.Stat: the score Z test statistic
- PVE: an approximation of the proportion of phenotype variance explained

Unused RVTESTS fields:
- INFORMATIVE_ALT_AC: The number of alternative alleles in the analyzed samples.
- CALL_RATE: The fraction of non-missing alleles.
- HWE_PVALUE: Hardy-Weinberg equilibrium.
- N_REF: Number of samples carrying homozygous reference alleles.
- N_HET: Number of samples carrying heterozygous alleles.
- N_ALT: Number of samples carrying homozygous alternative alleles.
- U_STAT: Score statistic.
