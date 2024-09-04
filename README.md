
<!-- README.md is generated from README.Rmd. Please edit that file -->

# trastools

<!-- badges: start -->
<!-- badges: end -->

The goal of trastools is to provide simple functions that help parsing
info from transcriptome assembly files, such as GTF and tracking files
from gffCompare.

## Installation

You can install the development version of trastools from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("llorenzi90/trastools")
```

## Example

``` r
library(trastools)
tracking <- structure(list(
V1 = c("TCONS_00000001", "TCONS_00000002", "TCONS_00000003",
"TCONS_00000004", "TCONS_00000005", "TCONS_00000006"),
V2 = c("XLOC_000001",
       "XLOC_000003",
       "XLOC_000004",
       "XLOC_000005",
       "XLOC_000009",
       "XLOC_000009"
),
V3 = c("XLOC_002036|TCONS_00004598", "-", "-", "XLOC_002044|TCONS_00004614",
          "XLOC_000074|TCONS_00000075", "XLOC_000074|TCONS_00000080"),
V4 = c("i", "u", "r", "x", "=", "="),
V5 = c("q1:STRG.1|STRG.1.1|1|0.091235|0.267011|6.848781|410",
       "q1:STRG.5|STRG.5.1|1|0.066681|0.195150|5.005556|360",
       "q1:STRG.10|STRG.10.1|1|0.084226|0.246497|6.322581|248",
       "q1:STRG.13|STRG.13.1|1|0.079929|0.233920|6.000000|251",
       "q1:STRG.32|STRG.32.1|9|8.963122|26.231630|672.835327|2474",
       "q1:STRG.32|STRG.32.2|8|0.459790|1.345631|34.515125|2372"
       ),
V6 = c("-", "-", "-", "-",
"q2:STRG.11|STRG.11.1|9|8.893389|29.440718|746.519531|2430",
"q2:STRG.11|STRG.11.2|8|0.305824|1.012402|25.671185|2328"
),
V7 = c("-", "-", "-", "-",
"q3:STRG.23|STRG.23.2|9|8.601234|23.619555|632.080566|2429",
"q3:STRG.23|STRG.23.1|8|0.289132|0.793974|21.247471|2594"
)),
row.names = c(NA, 6L), class = "data.frame")

get_sample_occurrance_per_gene_from_tracking(tracking)
#>     V1    V2    V3     gene_id
#> 1 TRUE FALSE FALSE XLOC_000001
#> 2 TRUE FALSE FALSE XLOC_000003
#> 3 TRUE FALSE FALSE XLOC_000004
#> 4 TRUE FALSE FALSE XLOC_000005
#> 5 TRUE  TRUE  TRUE XLOC_000009
```
