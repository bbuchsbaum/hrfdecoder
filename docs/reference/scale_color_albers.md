# Albers Color Scale for ggplot2

Discrete or continuous color scale matching the Albers minimalist theme.

## Usage

``` r
scale_color_albers(..., discrete = TRUE)
```

## Arguments

- ...:

  Additional arguments passed to scale functions

- discrete:

  Logical; if TRUE uses Okabe-Ito palette, if FALSE uses viridis

## Value

A ggplot2 scale object

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
  geom_point() +
  scale_color_albers()
} # }
```
