#' Okabe-Ito Color Palette for Albers Theme
#'
#' Returns the colorblind-safe Okabe-Ito palette for use in plots.
#'
#' @return A character vector of 8 hex colors
#' @keywords internal
albers_okabe_ito <- function() {
  c("#000000", "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
}

#' Albers Color Scale for ggplot2
#'
#' Discrete or continuous color scale matching the Albers minimalist theme.
#'
#' @param ... Additional arguments passed to scale functions
#' @param discrete Logical; if TRUE uses Okabe-Ito palette, if FALSE uses viridis
#' @return A ggplot2 scale object
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
#'   geom_point() +
#'   scale_color_albers()
#' }
scale_color_albers <- function(..., discrete = TRUE) {
  if (discrete) {
    ggplot2::scale_color_manual(values = albers_okabe_ito(), ...)
  } else if (requireNamespace("viridisLite", quietly = TRUE)) {
    ggplot2::scale_color_gradientn(colours = viridisLite::viridis(256), ...)
  } else {
    ggplot2::scale_color_gradient(low = "#cbd5e1", high = "#1f6feb", ...)
  }
}

#' Albers Fill Scale for ggplot2
#'
#' Discrete or continuous fill scale matching the Albers minimalist theme.
#'
#' @param ... Additional arguments passed to scale functions
#' @param discrete Logical; if TRUE uses Okabe-Ito palette, if FALSE uses viridis
#' @return A ggplot2 scale object
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(x = factor(cyl), fill = factor(cyl))) +
#'   geom_bar() +
#'   scale_fill_albers()
#' }
scale_fill_albers <- function(..., discrete = TRUE) {
  if (discrete) {
    ggplot2::scale_fill_manual(values = albers_okabe_ito(), ...)
  } else if (requireNamespace("viridisLite", quietly = TRUE)) {
    ggplot2::scale_fill_gradientn(colours = viridisLite::viridis(256), ...)
  } else {
    ggplot2::scale_fill_gradient(low = "#cbd5e1", high = "#1f6feb", ...)
  }
}

#' Albers Minimalist ggplot2 Theme
#'
#' A clean, minimal ggplot2 theme with subtle grid lines, top legend placement,
#' and accessible typography. Designed to match the Albers visual system used
#' in package vignettes.
#'
#' @param base_size Base font size in points (default: 11)
#' @param base_family Base font family (default: "system-ui")
#' @return A ggplot2 theme object
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#' theme_set(theme_albers())
#'
#' ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
#'   geom_point(size = 2.2) +
#'   scale_color_albers() +
#'   labs(
#'     title = "Fuel efficiency vs. weight",
#'     subtitle = "Example with Albers theme",
#'     x = "Weight (1000 lbs)",
#'     y = "MPG"
#'   )
#' }
theme_albers <- function(base_size = 11, base_family = "system-ui") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "#e5e7eb", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 6)),
      plot.title = ggplot2::element_text(face = "bold", margin = ggplot2::margin(b = 8)),
      plot.subtitle = ggplot2::element_text(margin = ggplot2::margin(b = 10), color = "#374151"),
      plot.caption = ggplot2::element_text(size = ggplot2::rel(0.9), color = "#6b7280",
                                           margin = ggplot2::margin(t = 10)),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold")
    )
}
