# Run from the repository root with: make -C notes
# Generated files stay in the ignored agent/ workspace.
LATEXMK = latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error \
	-emulate-aux-dir -auxdir="../agent/tmp/latex-notes" \
	-outdir="../agent/outputs/notes"
PYTHON ?= $(firstword $(wildcard $(HOME)/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3) python3)

.PHONY: pdf watch clean plot-labels
.DEFAULT_GOAL := pdf

pdf: plot-labels
	$(LATEXMK) notes.tex

plot-labels:
	@if test -f ../agent/scripts/prepare_notes_plot_labels.py; then \
		$(PYTHON) ../agent/scripts/prepare_notes_plot_labels.py; \
	else \
		mkdir -p ../agent/outputs/notes/plots; \
		cp ../output/pdf/multiqt_tuning.pdf ../output/pdf/multiqt_method_forecasts.pdf \
			../output/pdf/multiqt_weight_timelines.pdf ../output/pdf/multiqt_peak_reasoning_tradeoff.pdf \
			../output/pdf/multiqt_evaluation_tree.pdf ../output/pdf/multiqt_soft_switch_timelines.pdf \
			../agent/outputs/notes/plots/; \
		for group in multiqt_weight_timelines multiqt_method_forecasts; do \
			mkdir -p "../agent/outputs/notes/plots/$$group"; \
			cp "../output/pdf/$$group/flu_horizon2_delay0.pdf" \
				"../agent/outputs/notes/plots/$$group/"; \
		done; \
	fi

# Rebuild after each save; open the PDF separately. Stop with Ctrl-C.
watch: plot-labels
	$(LATEXMK) -pvc -view=none notes.tex

# Remove temporary compilation files; keep the source and finished PDF.
clean:
	$(LATEXMK) -c notes.tex
