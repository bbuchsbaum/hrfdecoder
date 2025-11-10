# Albers Minimalist ggplot2 Theme

A clean, minimal ggplot2 theme with subtle grid lines, top legend
placement, and accessible typography. Designed to match the Albers
visual system used in package vignettes.

## Usage

``` r
theme_albers(base_size = 11, base_family = "system-ui")
```

## Arguments

- base_size:

  Base font size in points (default: 11)

- base_family:

  Base font family (default: "system-ui")

## Value

A ggplot2 theme object

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
theme_set(theme_albers())

ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
  geom_point(size = 2.2) +
  scale_color_albers() +
  labs(
    title = "Fuel efficiency vs. weight",
    subtitle = "Example with Albers theme",
    x = "Weight (1000 lbs)",
    y = "MPG"
  )
} # }
```
