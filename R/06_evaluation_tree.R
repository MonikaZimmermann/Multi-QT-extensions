# Purpose: an overview of where the recalibration methods gain against the original
# forecasts. One tree per disease (covid and flu stay separate): the root pools all
# horizons, feedback designs and locations, the next levels split by feedback design and
# then by horizon. Every node shows three bar charts, overall, up to peak and after peak:
# the loss reduction against the original forecasts, each with the coverage error of the
# same methods below it. First page per disease: all locations except the US pooled
# (pooled_rows()); next pages: the locations of multiqt_settings() (US, California, Texas,
# New York, Pennsylvania), one location each. Each page draws its own best methods.
# Run from the project root after R/04_methods.R, which scores every location, horizon and
# feedback design of both diseases: Rscript R/06_evaluation_tree.R
# Order of the workflow: 00 -> 01 -> 05 -> 02 / 04 -> 06.
source("R/03_functions.R")

# 1. SETTINGS -----------------------------------------------------------------
settings <- multiqt_settings()
prepared <- readRDS("output/data/peak_aligned.rds")
scores <- read.csv("output/results/peak_reasoning/summary_scores_all_locations.csv",
  colClasses = c(location = "character"))
# Scores and prepared data of another feedback protocol (R/03_functions.R) are not combined.
check_feedback_protocol(scores$feedback_protocol, "summary_scores_all_locations.csv",
  "Rscript R/04_methods.R")
check_feedback_protocol(prepared$settings$feedback_protocol$id, "output/data/peak_aligned.rds",
  "Rscript R/01_data_overview.R")
diseases <- c("covid", "flu")
missing <- setdiff(diseases, scores$disease)
if (length(missing)) message("No scores for ", paste(missing, collapse = ", "),
  ": run R/04_methods.R (score_diseases).")
diseases <- setdiff(diseases, missing)
# Split order of the tree, top to bottom; swap the two to split by horizon first.
splits <- c("feedback", "horizon")
split_titles <- c(feedback = "feedback design", horizon = "horizon")
phases <- c("overall", "before", "after")
min_cases <- 5L                    # Nodes with fewer cases show no share of improved cases.
# Scored on the tracker's scale, as the tuning and the trade-off pages of R/04.
loss_column <- if (settings$transform == "log1p") "log_quantile_loss" else "quantile_loss"
loss_label <- if (settings$transform == "log1p") "loss of log(counts + 1)" else
  paste("loss in", count_scale_label(settings$scale_method))
# Methods drawn (the csv has every method): on each page the n_tree_methods with the
# largest loss reduction in the window selection_phase (after peak) at the root of that page
# (all its feedback designs, horizons and locations; the same original loss for every
# method, so the smallest pooled loss), among the methods of 04's plots (method_groups in
# R/03) except the target oracle (hindsight); best first, so the pages of a disease can draw
# different methods. Chosen on the evaluation period itself: a ranking, not a choice of
# method.
n_tree_methods <- 10L
selection_phase <- "after"
tree_candidates <- setdiff(unique(unlist(method_groups)), c("original", "oracle_target"))
# Colours as in 04's trade-off plots (tradeoff_colours in R/03), so a method keeps its
# colour on every page. A variant shares its method's shade: pale with an outline in the
# colour = common weights for all levels; black outline = + outer MultiQT. marker: the
# colour of the "<" / ">" on a bar cut at the axis, black or white, whichever contrasts more
# with the fill as drawn (WCAG luminance; a pale fill is its colour over the white panel).
bar_style <- function(methods) {
  colour <- tradeoff_colours[methods]
  common <- grepl("_common", methods)
  fill <- ifelse(common, adjustcolor(colour, alpha.f = .35), colour)
  rgba <- col2rgb(fill, alpha = TRUE) / 255
  shown <- t(rgba[1:3, , drop = FALSE]) * rgba[4, ] + 1 - rgba[4, ]
  luminance <- (ifelse(shown <= .04045, shown / 12.92, ((shown + .055) / 1.055)^2.4) %*%
    c(.2126, .7152, .0722))[, 1]
  list(fill = fill, border = ifelse(grepl("_outer$", methods), "black",
    ifelse(common, colour, NA)), marker = ifelse(luminance > .179, "black", "white"))
}
# The key's note on the colours and on the variants drawn, one line each, so it fits.
key_note <- function(methods) c("Colours as in the trade-off plots of R/04.",
  if (any(grepl("_common", methods))) "Pale: common weights for all levels.",
  if (any(grepl("_outer$", methods))) "Black outline: + outer MultiQT.")
window_titles <- c(overall = "Overall",
  before = paste0("Up to peak (weeks -", settings$weeks_before_peak - 1, " to 0)"),
  after = paste0("After peak (weeks 1 to ", settings$weeks_after_peak, ")"))
result_dir <- "output/results/peak_reasoning"
figure_dir <- "output/pdf"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# 2. LOSS REDUCTION AND COVERAGE ERROR PER NODE -------------------------------
# A case is one location, horizon and feedback design, scored over the evaluation seasons
# together. Feedback designs: "delay d" (learn once from the report d weeks after the first
# one) or "re-run" (re-run every week on the reports known then).
scores <- scores[scores$season == "all" & scores$phase %in% phases &
  is.finite(scores[[loss_column]]), ]
scores$feedback <- ifelse(scores$feedback_revisions, ifelse(scores$feedback_delay > 0,
  paste("re-run, delay", scores$feedback_delay), "re-run"), paste("delay", scores$feedback_delay))
feedback_levels <- unique(scores$feedback[order(scores$feedback_revisions, scores$feedback_delay)])
scores$loss <- scores[[loss_column]]
case_columns <- c("disease", "location", "horizon", "feedback", "phase")
original <- scores[scores$method == "original", c(case_columns, "loss", "coverage_error")]
names(original)[-seq_along(case_columns)] <- c("original_loss", "original_coverage_error")
scores <- merge(scores[scores$method != "original", ], original, by = case_columns)
split_values <- list(feedback = feedback_levels, horizon = sort(unique(scores$horizon)))
split_names <- function(split, value) if (split == "horizon") paste("horizon", value) else value

# Reduction in %: 100 * (1 - pooled loss / pooled original loss), each case's mean loss
# weighted by its number of weeks; share_improved: the share of cases in which the method's
# loss is below the original's. Coverage error (per case the mean over the levels of
# |coverage - alpha| over the weeks of the window) is pooled the same way, weighted by
# the weeks, for the method and for the original.
node_summary <- function(rows) do.call(rbind, lapply(split(rows, list(rows$method,
  rows$phase), drop = TRUE), function(s) data.frame(method = s$method[1], phase = s$phase[1],
    cases = nrow(s), weeks = sum(s$n_weeks),
    reduction = 100 * (1 - sum(s$n_weeks * s$loss) / sum(s$n_weeks * s$original_loss)),
    share_improved = mean(s$loss < s$original_loss),
    coverage_error = weighted.mean(s$coverage_error, s$n_weeks),
    original_coverage_error = weighted.mean(s$original_coverage_error, s$n_weeks))))
# The nodes of one tree, root first: level, the values of the splits above (NA: pooled),
# its row on the page (leaves fill the rows in order; a parent sits beside the middle of
# its children) and the summary of its cases.
tree_nodes <- function(rows) {
  n <- lengths(split_values[splits])
  grid <- expand.grid(rev(lapply(split_values[splits], seq_along)))[, rev(seq_along(splits))]
  nodes <- list()
  for (level in 0:length(splits)) {
    keys <- unique(grid[, seq_len(level), drop = FALSE])
    for (i in seq_len(if (level == 0) 1L else nrow(keys))) {
      leaves <- if (level == 0) seq_len(nrow(grid)) else which(apply(grid[, seq_len(level),
        drop = FALSE], 1, function(row) all(row == unlist(keys[i, ]))))
      key <- setNames(rep(NA_character_, length(splits)), splits)
      inside <- rep(TRUE, nrow(rows))
      for (j in seq_len(level)) {
        key[j] <- as.character(split_values[[splits[j]]][keys[i, j]])
        inside <- inside & as.character(rows[[splits[j]]]) == key[j]
      }
      nodes[[length(nodes) + 1L]] <- list(level = level, key = key,
        row = leaves[ceiling(length(leaves) / 2)], parent_row = NA,
        summary = if (any(inside)) node_summary(rows[inside, ]))
    }
  }
  # Each node's parent: the node one level up with the same split values above.
  for (k in seq_along(nodes)) if (nodes[[k]]$level > 0) {
    level <- nodes[[k]]$level
    parent <- Find(function(node) node$level == level - 1L &&
      identical(node$key[seq_len(level - 1L)], nodes[[k]]$key[seq_len(level - 1L)]), nodes)
    nodes[[k]]$parent_row <- parent$row
  }
  nodes
}
# Pages: NA = the pooled locations, then each location of settings$location_ids.
pages <- do.call(rbind, lapply(diseases, function(disease)
  data.frame(disease, location = c(NA, settings$location_ids))))
pages <- pages[is.na(pages$location) |                # Locations without scores: no page.
  paste(pages$disease, pages$location) %in% paste(scores$disease, scores$location), ]
page_rows <- lapply(seq_len(nrow(pages)), function(p) {
  rows <- scores[scores$disease == pages$disease[p], ]
  if (is.na(pages$location[p])) pooled_rows(rows, settings) else
    rows[rows$location == pages$location[p], ]
})
trees <- lapply(page_rows, tree_nodes)
# The methods of each page (see n_tree_methods), with their labels.
tree_methods <- lapply(trees, function(nodes) {
  root <- nodes[[1]]$summary
  root <- root[root$phase == selection_phase & root$method %in% tree_candidates, ]
  top <- head(root$method[order(-root$reduction)], n_tree_methods)
  setNames(method_short_labels[top], top)
})
write.csv(do.call(rbind, lapply(seq_len(nrow(pages)), function(p)
  do.call(rbind, lapply(trees[[p]], function(node) if (!is.null(node$summary))
    data.frame(disease = pages$disease[p], location = ifelse(is.na(pages$location[p]), "all",
      pages$location[p]), level = node$level, as.list(node$key), node$summary,
      feedback_protocol = feedback_protocol()$id))))),
  file.path(result_dir, "evaluation_tree.csv"), row.names = FALSE)
# Footnote on the drawn methods with identical forecasts. In each case R/04 names, for a
# method whose forecasts equal an earlier one's (method_names order), the first such method
# (duplicate_of); two methods are identical in a case if they lead to the same first method.
# Each pair identical in all of the page's cases of a feedback design, with those designs
# (pairs identical only in some cases are left out: at horizon -1, without a target signal,
# Decline MultiQT gives the original forecast and score_gated that of score_learned, which
# the second line says once).
twin_note <- function(rows, methods) {
  rows <- rows[rows$phase == "overall" & rows$method %in% names(methods), ]
  rows$first <- ifelse(is.na(rows$duplicate_of), rows$method, rows$duplicate_of)
  twins <- merge(rows, rows, by = c("location", "horizon", "feedback", "first"))
  twins <- twins[match(twins$method.x, method_names) > match(twins$method.y, method_names), ]
  cases <- tapply(paste(rows$location, rows$horizon), rows$feedback,
    function(case) length(unique(case)))
  pairs <- split(twins$feedback, paste(methods[twins$method.x], "=", methods[twins$method.y]))
  pairs <- Filter(length, lapply(pairs, function(feedback) {
    n <- table(factor(feedback, feedback_levels))
    names(n)[n > 0 & n == cases[names(n)]]
  }))
  context <- intersect(names(methods), c("decline", "score_gated"))
  c(if (length(pairs)) paste("Identical forecasts in every case of a design:",
    paste(paste0(names(pairs), " (", vapply(pairs, paste, "", collapse = ", "), ")"),
      collapse = "; ")),
    if (length(context)) paste("At horizon -1 (no target signal) Decline MultiQT gives the",
      "original forecast and the learned score mixture does not back off."))
}

# 3. PLOT ---------------------------------------------------------------------
# Columns: the three windows of each tree level, with a gap for the connectors between
# levels; rows: the leaves, each a loss panel with the coverage error (CE) panel of its
# window below. The method key fills the free space above the root. All loss panels of a
# page share one axis: up to 1.15 times the 95th percentile of the absolute reductions; a
# bar beyond it ends in "<" or ">". All CE panels of a page share another, from 0.
panel_width <- 2.05
gap_width <- .55
# Inches per leaf: the loss panel and the CE panel below it. The loss panel's margins are
# 1.2 lines (.16 inches) larger, for the titles above it, so the bars of both panels get
# the same height.
loss_height <- .98
ce_height <- .82
n_levels <- length(splits) + 1L
n_rows <- prod(lengths(split_values[splits]))
column_of <- function(level, w) level * (length(phases) + 1L) + w
widths <- rep(c(rep(panel_width, length(phases)), gap_width), n_levels)[-n_levels *
  (length(phases) + 1L)]
# methods: named labels of the methods drawn, top to bottom; style: bar_style(methods).
draw_panel <- function(s, xlim, show_share, methods, style) {
  values <- s$reduction[match(names(methods), s$method)]    # NULL s: no weeks.
  shares <- s$share_improved[match(names(methods), s$method)]
  y <- seq_along(methods)
  par(mar = c(1.25, .5, 2, if (show_share) 2.2 else .8), mgp = c(1, .25, 0), tcl = -.2)
  plot(NA, xlim = xlim, ylim = c(length(y) + .6, .4), xaxs = "i", yaxs = "i", xaxt = "n",
    yaxt = "n", xlab = "", ylab = "")
  ticks <- pretty(xlim, n = 5)
  abline(v = ticks, col = "gray90")
  axis(1, at = ticks, labels = paste0(ticks, "%"), cex.axis = .6)
  drawn <- is.finite(values)
  if (!any(drawn)) {
    text(mean(xlim), mean(range(y)), "No weeks in this window", cex = .7, col = "gray40")
    return(invisible())
  }
  shown <- pmin(pmax(values, xlim[1]), xlim[2])
  rect(0, y - .38, shown, y + .38, col = style$fill, border = style$border, lwd = .7)
  # One call per side: text() uses only the first two values of adj.
  for (right in c(FALSE, TRUE)) {
    cut <- which(drawn & (if (right) values > xlim[2] else values < xlim[1]))
    if (length(cut)) text(shown[cut], y[cut], if (right) ">" else "<", adj = c(right, .5),
      cex = .7, font = 2, col = style$marker[cut])
  }
  abline(v = 0, col = "gray30")
  box(col = "gray60")
  if (show_share) mtext(ifelse(drawn, sprintf("%.0f%%", 100 * shares), ""), side = 4, at = y,
    las = 1, line = .15, cex = .38, col = "gray25")
}
# The CE panel below a loss panel: the same methods and bars, from 0; dashed: the
# original's CE pooled over all the node's cases. Its right margin is the loss panel's, so
# the two panels line up.
draw_ce_panel <- function(s, xlim, show_share, methods, style) {
  values <- s$coverage_error[match(names(methods), s$method)]
  y <- seq_along(methods)
  par(mar = c(1.2, .5, .85, if (show_share) 2.2 else .8), mgp = c(1, .25, 0), tcl = -.2)
  plot(NA, xlim = xlim, ylim = c(length(y) + .6, .4), xaxs = "i", yaxs = "i", xaxt = "n",
    yaxt = "n", xlab = "", ylab = "")
  ticks <- pretty(xlim, n = 4)
  abline(v = ticks, col = "gray90")
  axis(1, at = ticks, cex.axis = .6)
  mtext("Coverage error", side = 3, line = .05, adj = 0, cex = .45, col = "gray35")
  if (any(is.finite(values))) {
    rect(0, y - .38, values, y + .38, col = style$fill, border = style$border, lwd = .7)
    abline(v = s$original_coverage_error[which.max(s$cases)], lty = 2, lwd = .9)
  }
  box(col = "gray60")
}
page_title <- function(disease, location) {
  where <- if (is.na(location)) paste("all", pooled_label(length(unique(pooled_rows(
    scores[scores$disease == disease, ], settings)$location)), settings), "pooled") else {
    name <- unique(scores$location_name[scores$location == location])[1]
    paste0(name, " (", location, ")")
  }
  paste(disease, "|", where, "| evaluation:", prepared$settings$evaluation_label, "|",
    loss_label)
}
# methods: the page's methods (tree_methods); twins: its twin_note(), NULL without twins.
draw_tree_page <- function(nodes, methods, twins, disease, location) {
  style <- bar_style(names(methods))
  values <- unlist(lapply(nodes, function(node) node$summary$reduction[
    node$summary$method %in% names(methods)]))
  limit <- min(max(abs(values), na.rm = TRUE), 1.15 * quantile(abs(values), .95, na.rm = TRUE))
  xlim <- c(-1, 1) * limit
  if (min(values, na.rm = TRUE) > -limit / 4) xlim[1] <- -limit / 4   # Room for gains first.
  if (max(values, na.rm = TRUE) < limit / 4) xlim[2] <- limit / 4
  # CE axis: 0 to the largest CE drawn, of the methods or of the original.
  ce_values <- unlist(lapply(nodes, function(node) c(node$summary$coverage_error[
    node$summary$method %in% names(methods)], node$summary$original_coverage_error)))
  ce_xlim <- c(0, 1.04 * max(ce_values, na.rm = TRUE))
  # Panel numbers per node and window, the loss panel in a leaf's first layout row and its
  # CE panel in the second; 0: empty.
  index <- matrix(0L, 2L * n_rows, length(widths))
  for (k in seq_along(nodes)) for (w in seq_along(phases))
    index[2L * nodes[[k]]$row - 1:0, column_of(nodes[[k]]$level, w)] <-
      (k - 1L) * 2L * length(phases) + 2L * w - 1:0
  key_rows <- seq_len(2L * max(1L, nodes[[1]]$row - 2L))
  index[key_rows, column_of(0L, seq_along(phases))] <- 2L * length(nodes) * length(phases) + 1L
  layout(index, widths = widths, heights = rep(c(loss_height, ce_height), n_rows))
  par(oma = c(5.25, 1, 3.3, 1))
  edges <- list()
  for (k in seq_along(nodes)) {
    node <- nodes[[k]]
    for (w in seq_along(phases)) {
      s <- node$summary[node$summary$phase == phases[w], ]
      show_share <- !is.null(s) && nrow(s) > 0 && max(s$cases) >= min_cases
      draw_panel(s, xlim, show_share, methods, style)
      mtext(window_titles[phases[w]], side = 3, line = .15, cex = .5, col = "gray20")
      if (w == 2) {
        parts <- ifelse(is.na(node$key), paste0("all ", split_titles[splits], "s"),
          mapply(split_names, splits, node$key))
        cases <- if (is.null(node$summary)) 0 else max(node$summary$cases)
        mtext(paste0(paste(parts, collapse = " | "), " | ", cases,
          if (cases == 1) " case" else " cases"), side = 3, line = .95, cex = .62, font = 2)
      }
      # Connectors meet a node between its loss and CE panels.
      if (w == 1) edges[[k]] <- c(left = grconvertX(0, "nfc", "ndc"),
        y = grconvertY(0, "nfc", "ndc"))
      if (w == length(phases)) edges[[k]]["right"] <- grconvertX(1, "nfc", "ndc")
      draw_ce_panel(s, ce_xlim, show_share, methods, style)
    }
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  key <- legend("center", unname(methods), fill = style$fill, border = style$border,
    bty = "n", cex = 1.2, y.intersp = 1.3,
    title = paste("Bars, top to bottom: the", length(methods), "best on this page",
      paste0("(quantile loss ", tolower(sub(" \\(.*", "", window_titles[[selection_phase]])),
        ", all designs and horizons)")), title.font = 2)
  text(key$rect$left, key$rect$top - key$rect$h, paste(key_note(names(methods)),
    collapse = "\n"), adj = c(0, 1), cex = .95, col = "gray25")
  # Connectors from each node to its children, on top of the whole page.
  par(oma = c(0, 0, 0, 0))
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  rows <- vapply(nodes, `[[`, numeric(1), "row")
  levels <- vapply(nodes, `[[`, numeric(1), "level")
  for (k in seq_along(nodes)) if (nodes[[k]]$level > 0) {
    parent <- edges[[which(levels == nodes[[k]]$level - 1L & rows == nodes[[k]]$parent_row)]]
    child <- edges[[k]]
    middle <- (parent["right"] + child["left"]) / 2
    segments(c(parent["right"], middle, middle), c(parent["y"], parent["y"], child["y"]),
      c(middle, middle, child["left"]), c(parent["y"], child["y"], child["y"]),
      col = "gray45", lwd = 1.2)
  }
  mtext(page_title(disease, location), side = 3, line = -2, font = 2, cex = 1.35)
  notes <- c(paste("Loss bar: 100 * (1 - pooled loss / pooled loss of the original), pooled over",
    "the node's cases (location x horizon x feedback design), each case's mean loss weighted by",
    "its number of weeks; to the right: positive = smaller loss than the original."),
    paste0("Numbers right of each loss panel: share of the node's cases in which the method's ",
      "loss is below the original's (nodes with at least ", min_cases, " cases). Up to peak / ",
      "after peak: the weeks around every local peak of the latest counts."),
    paste("Below each loss panel: coverage error, the mean over the quantile levels of",
      "|coverage - alpha| over the weeks of the window, averaged over the cases of the node",
      "weighted by their weeks; dashed: the original forecasts."),
    twins)
  mtext(notes, side = 1, line = -4.25 + 1.05 * (seq_along(notes) - 1), cex = .8)
}
# Height: the leaves plus the outer margins (oma, 8.55 lines of .132 inches).
cairo_pdf(file.path(figure_dir, "multiqt_evaluation_tree.pdf"), width = sum(widths) + .5,
  height = n_rows * (loss_height + ce_height) + 1.15, onefile = TRUE,
  family = "DejaVu Sans")
for (p in seq_len(nrow(pages))) draw_tree_page(trees[[p]], tree_methods[[p]],
  twin_note(page_rows[[p]], tree_methods[[p]]), pages$disease[p], pages$location[p])
dev.off()
