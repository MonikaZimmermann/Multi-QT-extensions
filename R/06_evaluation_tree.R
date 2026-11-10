# Purpose: an overview of where the recalibration methods gain against the original
# forecasts. One tree per disease (covid and flu stay separate): the root pools all
# horizons, feedback designs and locations, the next levels split by feedback design and
# then by horizon. Every node shows three bar charts, overall, up to peak and after peak:
# the loss reduction against the original forecasts. First page per disease: all
# locations except the US pooled (pooled_rows()); next pages: the locations of
# multiqt_settings() (US, California, Texas, New York, Pennsylvania), one location each.
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
# Methods drawn, one or two per family (the csv has every method), in this order; Delphi's
# production trackers (covid: REF-op, flu: sqrt constant; Delphi's constants) as benchmarks.
tree_methods <- c(fixed = "Fixed MultiQT", ref_op = "Delphi REF-op (covid prod)",
  sqrt_constant = "Delphi sqrt 0.018 (flu prod)", studentized = "Studentized median",
  studentized_levels = "Studentized full base (level scales)",
  weighted_common = "Rate grid: learned, common weights",
  shared_common = "Rate grid: learned + sharing, common weights",
  restart_weighted = "Weekly starts: learned", diverse_uniform = "Diverse experts: equal",
  diverse_soft_switch = "Diverse experts: soft switching")
method_colours <- c(fixed = "#0072B2", ref_op = "#B2182B", sqrt_constant = "#A6761D",
  studentized = "#D55E00",
  studentized_levels = "#009E73",
  weighted_common = "#56B4E9", shared_common = "#882255", restart_weighted = "#999933",
  diverse_uniform = "#E69F00", diverse_soft_switch = "#7B3294")
window_titles <- c(overall = "Overall",
  before = paste0("Up to peak (weeks -", settings$weeks_before_peak - 1, " to 0)"),
  after = paste0("After peak (weeks 1 to ", settings$weeks_after_peak, ")"))
result_dir <- "output/results/peak_reasoning"
figure_dir <- "output/pdf"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# 2. LOSS REDUCTION PER NODE --------------------------------------------------
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
original <- scores[scores$method == "original", c(case_columns, "loss")]
names(original)[names(original) == "loss"] <- "original_loss"
scores <- merge(scores[scores$method != "original", ], original, by = case_columns)
split_values <- list(feedback = feedback_levels, horizon = sort(unique(scores$horizon)))
split_names <- function(split, value) if (split == "horizon") paste("horizon", value) else value

# Reduction in %: 100 * (1 - pooled loss / pooled original loss), each case's mean loss
# weighted by its number of weeks; share_improved: the share of cases in which the method's
# loss is below the original's.
node_summary <- function(rows) do.call(rbind, lapply(split(rows, list(rows$method,
  rows$phase), drop = TRUE), function(s) data.frame(method = s$method[1], phase = s$phase[1],
    cases = nrow(s), weeks = sum(s$n_weeks),
    reduction = 100 * (1 - sum(s$n_weeks * s$loss) / sum(s$n_weeks * s$original_loss)),
    share_improved = mean(s$loss < s$original_loss))))
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
trees <- lapply(seq_len(nrow(pages)), function(p) {
  rows <- scores[scores$disease == pages$disease[p], ]
  rows <- if (is.na(pages$location[p])) pooled_rows(rows, settings) else
    rows[rows$location == pages$location[p], ]
  tree_nodes(rows)
})
write.csv(do.call(rbind, lapply(seq_len(nrow(pages)), function(p)
  do.call(rbind, lapply(trees[[p]], function(node) if (!is.null(node$summary))
    data.frame(disease = pages$disease[p], location = ifelse(is.na(pages$location[p]), "all",
      pages$location[p]), level = node$level, as.list(node$key), node$summary,
      feedback_protocol = feedback_protocol()$id))))),
  file.path(result_dir, "evaluation_tree.csv"), row.names = FALSE)

# 3. PLOT ---------------------------------------------------------------------
# Columns: the three windows of each tree level, with a gap for the connectors between
# levels; rows: the leaves. The method key fills the free space above the root. All panels of a page share one axis: up to 1.15 times the 95th
# percentile of the absolute reductions; a bar beyond it ends in "<" or ">".
panel_width <- 2.05
gap_width <- .55
row_height <- 1.3
n_levels <- length(splits) + 1L
n_rows <- prod(lengths(split_values[splits]))
column_of <- function(level, w) level * (length(phases) + 1L) + w
widths <- rep(c(rep(panel_width, length(phases)), gap_width), n_levels)[-n_levels *
  (length(phases) + 1L)]
draw_panel <- function(s, xlim, show_share) {
  values <- s$reduction[match(names(tree_methods), s$method)]    # NULL s: no weeks.
  shares <- s$share_improved[match(names(tree_methods), s$method)]
  y <- seq_along(tree_methods)
  par(mar = c(1.7, .5, 2.2, if (show_share) 2.2 else .8), mgp = c(1, .25, 0), tcl = -.2)
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
  rect(0, y - .38, shown, y + .38, col = method_colours[names(tree_methods)], border = NA)
  beyond <- drawn & (values < xlim[1] | values > xlim[2])
  if (any(beyond)) text(shown[beyond], y[beyond], ifelse(values[beyond] > 0, ">", "<"),
    adj = c(ifelse(values[beyond] > 0, 1, 0), .5), cex = .8, font = 2, col = "white")
  abline(v = 0, col = "gray30")
  box(col = "gray60")
  if (show_share) mtext(ifelse(drawn, sprintf("%.0f%%", 100 * shares), ""), side = 4, at = y,
    las = 1, line = .15, cex = .42, col = "gray25")
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
draw_tree_page <- function(nodes, disease, location) {
  values <- unlist(lapply(nodes, function(node) node$summary$reduction[
    node$summary$method %in% names(tree_methods)]))
  limit <- min(max(abs(values), na.rm = TRUE), 1.15 * quantile(abs(values), .95, na.rm = TRUE))
  xlim <- c(-1, 1) * limit
  if (min(values, na.rm = TRUE) > -limit / 4) xlim[1] <- -limit / 4   # Room for gains first.
  if (max(values, na.rm = TRUE) < limit / 4) xlim[2] <- limit / 4
  # Panel numbers per node and window; 0: empty.
  index <- matrix(0L, n_rows, length(widths))
  for (k in seq_along(nodes)) for (w in seq_along(phases))
    index[nodes[[k]]$row, column_of(nodes[[k]]$level, w)] <- (k - 1L) * length(phases) + w
  key_rows <- seq_len(max(1L, nodes[[1]]$row - 2L))
  index[key_rows, column_of(0L, seq_along(phases))] <- length(nodes) * length(phases) + 1L
  layout(index, widths = widths, heights = rep(row_height, n_rows))
  par(oma = c(5.25, 1, 6.2, 1))
  edges <- list()
  for (k in seq_along(nodes)) {
    node <- nodes[[k]]
    for (w in seq_along(phases)) {
      s <- node$summary[node$summary$phase == phases[w], ]
      draw_panel(s, xlim, show_share = !is.null(s) && nrow(s) > 0 && max(s$cases) >= min_cases)
      mtext(window_titles[phases[w]], side = 3, line = .15, cex = .5, col = "gray20")
      if (w == 2) {
        parts <- ifelse(is.na(node$key), paste0("all ", split_titles[splits], "s"),
          mapply(split_names, splits, node$key))
        cases <- if (is.null(node$summary)) 0 else max(node$summary$cases)
        mtext(paste0(paste(parts, collapse = " | "), " | ", cases,
          if (cases == 1) " case" else " cases"), side = 3, line = 1, cex = .62, font = 2)
      }
      if (w == 1) edges[[k]] <- c(left = grconvertX(0, "nfc", "ndc"),
        y = grconvertY(.5, "npc", "ndc"))
      if (w == length(phases)) edges[[k]]["right"] <- grconvertX(1, "nfc", "ndc")
    }
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", unname(tree_methods), fill = method_colours[names(tree_methods)],
    border = NA, bty = "n", cex = 1.2, y.intersp = 1.3,
    title = "Bars, top to bottom", title.font = 2)
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
  mtext(c(paste("Loss reduction against the original forecasts. Left: all feedback designs and",
    "horizons pooled; then split by", paste0(paste(split_titles[splits],
      collapse = ", then by "), ".")), paste0("Feedback design (protocol ",
    feedback_protocol()$id, "): delay d = learn once, d weeks after the first report, from ",
    "the report current then; re-run = re-run every week on the reports known then.")),
    side = 3, line = c(-3.7, -4.9), cex = .95)
  notes <- c(paste("Bar: 100 * (1 - pooled loss / pooled loss of the original), pooled over the",
    "node's cases (location x horizon x feedback design), each case's mean loss weighted by",
    "its number of weeks; to the right: positive = smaller loss than the original."),
    paste0("Numbers right of each panel: share of the node's cases in which the method's loss is below ",
      "the original's (nodes with at least ", min_cases, " cases). Up to peak / after peak: ",
      "the weeks around every local peak of the latest counts."),
    paste("Each disease and feedback design uses the settings tuned with it",
      "(output/results/tuning); every method: output/results/peak_reasoning/evaluation_tree.csv."),
    paste("Issue cutoffs: Wednesdays (actual generation dates where known). Every design is",
      "scored against the latest reports, version", format(prepared$settings$truth_version),
      "(not against the reports learned from)."))
  mtext(notes, side = 1, line = -4.25 + 1.05 * (seq_along(notes) - 1), cex = .8)
}
cairo_pdf(file.path(figure_dir, "multiqt_evaluation_tree.pdf"), width = sum(widths) + .5,
  height = n_rows * row_height + 2.3, onefile = TRUE,
  family = "DejaVu Sans")
for (p in seq_len(nrow(pages))) draw_tree_page(trees[[p]], pages$disease[p], pages$location[p])
dev.off()
