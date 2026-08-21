# Right whale annual detection probability in Cape Cod Bay

Table S2 of Ganley et al. (2019): the goodness-of-fit and average
detection probability of a separately fitted detection function for each
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

## This is `p`, not perception bias

`p` here is the average detection probability of the fitted detection
function — the same quantity
[`selection_table()`](https://camilleross.org/dsfit/reference/selection_table.md)
reports as `p`, with its standard error as `p_se`. It is **not** a
`g(0)` component. Ganley et al. say so directly: perception bias "was
not addressed directly in this study", because estimating it needs a
second observer team they did not have.

The distinction matters because the two are easy to conflate and
correcting for one while believing you have corrected for the other is
how these estimates go wrong by a factor. See
[`g0()`](https://camilleross.org/dsfit/reference/g0.md), which will name
any component you have not supplied.

## What it is good for

It is a long, real example of the thing
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
is built around: **detection probability is not a constant of the
survey**. Twenty years of the same programme, the same aircraft and the
same bay give `p` between 0.431 and 0.866 — a twofold range — with
standard errors spanning an order of magnitude, from 0.028 to 0.333.

The goodness-of-fit columns are worth reading next to it. Ganley et al.
took p \> 0.05 as adequate fit, and 2003 fails on Kolmogorov-Smirnov
(0.048) while passing Cramér-von Mises (0.103); it also carries much the
largest standard error on `p`. A model can top a ranking and still fail
a fit test, which is why
[`selection_table()`](https://camilleross.org/dsfit/reference/selection_table.md)
carries `cvm_p` alongside the AIC.

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

# Detection probability is not a constant of a survey programme
range(ganley_detection$p)
#> [1] 0.431 0.866

# The year that fails Kolmogorov-Smirnov also has the worst precision on p
ganley_detection[ganley_detection$ks_p < 0.05, ]
#> # A tibble: 1 × 5
#>    year cvm_p  ks_p     p  p_se
#>   <int> <dbl> <dbl> <dbl> <dbl>
#> 1  2003 0.103 0.048  0.76 0.333
```
