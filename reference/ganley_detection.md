# Right whale annual detection probability in Cape Cod Bay

Table S2 of Ganley et al. (2019): goodness-of-fit p-values and average
detection probability for a detection function fitted separately to each
year from 1998 to 2017, from aerial line-transect surveys of Cape Cod
Bay.

## Usage

``` r
ganley_detection
```

## Format

A tibble with 20 rows and 5 columns:

- year:

  Survey year, 1998 to 2017.

- cvm_p:

  Cramér-von Mises p-value for the year's detection function.

- ks_p:

  Kolmogorov-Smirnov p-value.

- p:

  Average detection probability, \\\hat{P}\_a\\.

- p_se:

  Standard error of `p`.

## Source

Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what
there is: estimating North Atlantic right whale *Eubalaena glacialis*
local abundance. *Endangered Species Research* 38:101-113.
[doi:10.3354/esr00938](https://doi.org/10.3354/esr00938) , Table S2.
Open access under CC-BY.

## This is `p`, not perception

`p` is the average detection probability of the fitted detection
function, the same quantity
[`selection_table()`](https://camilleross.org/dsfit/reference/selection_table.md)
reports as `p` (with standard error `p_se`). It is not a `g(0)`
component. Ganley et al. state that perception bias was not addressed in
their study, because estimating it needs a second observer team.
Correcting for one while believing you have corrected for the other
biases abundance by a factor.
[`g0()`](https://camilleross.org/dsfit/reference/g0.md) names any
standard component that is missing.

## What it is useful for

It is a long, real example that detection probability is not a constant
of a survey. Twenty years of the same programme, aircraft and bay give
`p` between 0.431 and 0.866, with standard errors from 0.028 to 0.333.

The goodness-of-fit columns show that a model can rank well and still
fail a fit test. Ganley et al. took p \> 0.05 as adequate fit. 2003
fails Kolmogorov-Smirnov (0.048) but passes Cramér-von Mises (0.103),
and has the largest standard error on `p`. This is why
[`selection_table()`](https://camilleross.org/dsfit/reference/selection_table.md)
reports `cvm_p` alongside AIC.

## See also

[`selection_table()`](https://camilleross.org/dsfit/reference/selection_table.md),
[ganley_availability](https://camilleross.org/dsfit/reference/ganley_availability.md)

## Examples

``` r
ganley_detection
#> # A tibble: 20 × 5
#>     year cvm_p  ks_p     p  p_se
#>    <int> <dbl> <dbl> <dbl> <dbl>
#>  1  1998 0.751 0.564 0.663 0.156
#>  2  1999 0.521 0.426 0.601 0.094
#>  3  2000 0.531 0.248 0.75  0.041
#>  4  2001 0.867 0.802 0.64  0.07 
#>  5  2002 0.97  0.965 0.708 0.104
#>  6  2003 0.103 0.048 0.76  0.333
#>  7  2004 0.969 0.965 0.655 0.142
#>  8  2005 0.985 0.995 0.866 0.059
#>  9  2006 0.431 0.392 0.703 0.163
#> 10  2007 0.866 0.777 0.712 0.065
#> 11  2008 0.97  0.977 0.595 0.05 
#> 12  2009 0.581 0.507 0.591 0.038
#> 13  2010 0.668 0.785 0.599 0.046
#> 14  2011 0.847 0.737 0.604 0.043
#> 15  2012 0.991 0.982 0.611 0.062
#> 16  2013 0.98  0.993 0.576 0.03 
#> 17  2014 0.93  0.858 0.507 0.035
#> 18  2015 0.583 0.815 0.541 0.028
#> 19  2016 0.115 0.065 0.431 0.038
#> 20  2017 0.341 0.129 0.518 0.028

# Detection probability varies from year to year
range(ganley_detection$p)
#> [1] 0.431 0.866

# The year that fails Kolmogorov-Smirnov also has the least precise p
ganley_detection[ganley_detection$ks_p < 0.05, ]
#> # A tibble: 1 × 5
#>    year cvm_p  ks_p     p  p_se
#>   <int> <dbl> <dbl> <dbl> <dbl>
#> 1  2003 0.103 0.048  0.76 0.333
```
