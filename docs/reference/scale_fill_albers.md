# Albers Fill Scale for ggplot2

Discrete or continuous fill scale matching the Albers minimalist theme.

## Usage

``` r
scale_fill_albers(..., discrete = TRUE)
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
ggplot(mtcars, aes(x = factor(cyl), fill = factor(cyl))) +
  geom_bar() +
  scale_fill_albers()
} # }
```
