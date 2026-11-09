#' Plot the profile of motif of specific peak
#'
#' @param df motif information data.frame.
#' @param legend_lab legend lab.
#' @param x_lab x axis label.
#' @param y_lab y axis label.
#' @param interactive produce interactive fig or not.
#' @param width_svg width_svg
#' @param height_svg height_svg
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 labs
#' @importFrom dplyr left_join
#' @importFrom dplyr summarise
#' @return ggplot object
#' @examples
#' require(BSgenome.Hsapiens.UCSC.hg38)
#' data(pwm_obj)
#' region_oi <- GRanges(
#'     seqnames = "chr22",
#'     ranges = IRanges(start = 10525891, end = 10525991)
#' )
#' motifMatrix <- getMotifMatrix(
#'     region = region_oi,
#'     pwm = pwm_obj[c(45, 120, 170)],
#'     ref_obj = BSgenome.Hsapiens.UCSC.hg38
#' )
#' plotMotifProf(motifMatrix)
#' @export
plotMotifProf <- function(df, legend_lab = "motif", y_lab = "motif score",
                          x_lab = NULL, interactive = FALSE, width_svg = 10, height_svg = 6) {
    range_x <- attr(df, "range")
    if (!is.null(df$region) && length(unique(df$chr)) > 1) {
        stop("`df` holds motif hits of regions on several chromosomes, which cannot share one axis; ",
             "plot one region at a time, e.g. plotMotifProf(getMotifMatrix(region[1], pwm, ref_obj)).")
    }
    if (!is.null(df$region) && length(unique(df$region)) > 1) {
        message("`df` holds motif hits of ", length(unique(df$region)), " regions; they are drawn on one axis ",
                "spanning ", range_x[1], "-", range_x[2], ".")
    }
    x_max <- range_x[2]
    x_min <- range_x[1]
    score_max <- round(max(abs(range(df$score))))
    if (is.null(x_lab)) {
        x_lab <- paste0("peak ", df[1, "chr"], ":", x_min, "~", x_max)
    }


    coordinate <- score <- motif <- chr <- NULL
    motif_start <- motif_end <- motif_chr <- motif_strand <- motif_score <- seg_id <- NULL

    if (interactive) {
        df_interactive <- df %>%
            group_by(motif, seg_id) %>%
            summarise(
                motif_chr = unique(chr),
                motif_start = min(coordinate),
                motif_strand = unique(strand),
                motif_end = max(coordinate),
                motif_score = unique(score),
                .groups = "drop"
            ) %>%
            left_join(df, by = c("motif", "seg_id"))

        rlang::check_installed("ggiraph", reason = "For interactive plot.")

        if (requireNamespace("ggiraph", quietly = TRUE)) {
            p <- ggplot(df_interactive, mapping = aes(
                x = coordinate, y = score, color = motif, data_id = motif,
                group = interaction(motif, score, seg_id),
                tooltip = paste0(
                    "Motif: ", motif, "\n",
                    "chr: ", motif_chr, "\n",
                    "strand: ", motif_strand, "\n",
                    "start: ", motif_start, "\n",
                    "end: ", motif_end, "\n",
                    "score: ", round(score, 3)
                )
            )) +
                ggiraph::geom_line_interactive(size = 1) +
                geom_hline(yintercept = 0, color = "black", linewidth = 1) +
                labs(color = legend_lab, y = y_lab, x = x_lab) +
                coord_cartesian(xlim = c(x_min, x_max), ylim = c(-score_max, score_max)) +
                scale_x_continuous(
                    breaks = round(as.numeric(stats::quantile(seq(x_min, x_max), c(0, 0.25, 0.5, 0.75, 1)))),
                    labels = round(as.numeric(stats::quantile(seq(x_min, x_max), c(0, 0.25, 0.5, 0.75, 1))))
                ) +
                scale_y_continuous(
                    breaks = round(c(score_max * (-1), score_max * (-0.5), 0, score_max * 0.5, score_max)),
                    labels = round(c(score_max * (-1), score_max * (-0.5), 0, score_max * 0.5, score_max))
                ) +
                theme_classic() +
                theme(
                    panel.grid.minor = element_blank(),
                    axis.line.x = element_line(colour = "black")
                )

            p <- ggiraph::girafe(
                ggobj = p, width_svg = width_svg, height_svg = height_svg,
                options = list(
                    ggiraph::opts_hover(css = "stroke: orange; stroke-width: 2px; cursor: pointer;"),
                    ggiraph::opts_hover_inv(css = "opacity: 0.5;"),
                    ggiraph::opts_zoom(min = .5, max = 5),
                    ggiraph::opts_tooltip(css = "background-color: white; border: 1px solid #333; padding: 8px; border-radius: 4px; font-size: 12px;")
                )
            )
        }
    } else {
        p <- ggplot(df, mapping = aes(x = coordinate, y = score, color = motif, group = interaction(motif, score, seg_id))) +
            geom_line() +
            geom_hline(yintercept = 0, color = "black", linewidth = 1) +
            labs(color = legend_lab, y = y_lab, x = x_lab) +
            coord_cartesian(xlim = c(x_min, x_max), ylim = c(-score_max, score_max)) +
            scale_x_continuous(
                breaks = round(as.numeric(stats::quantile(seq(x_min, x_max), c(0, 0.25, 0.5, 0.75, 1)))),
                labels = round(as.numeric(stats::quantile(seq(x_min, x_max), c(0, 0.25, 0.5, 0.75, 1))))
            ) +
            scale_y_continuous(
                breaks = round(c(score_max * (-1), score_max * (-0.5), 0, score_max * 0.5, score_max)),
                labels = round(c(score_max * (-1), score_max * (-0.5), 0, score_max * 0.5, score_max))
            ) +
            theme_classic() +
            theme(
                panel.grid.minor = element_blank(),
                axis.line.x = element_line(colour = "black")
            )
    }


    return(p)
}
