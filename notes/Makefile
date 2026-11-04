# Run from the repository root with: make -C notes
# Generated files stay in the ignored agent/ workspace.
LATEXMK = latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error \
	-emulate-aux-dir -auxdir="../agent/tmp/latex-notes" \
	-outdir="../agent/outputs/notes"

.PHONY: pdf watch clean
.DEFAULT_GOAL := pdf

pdf:
	$(LATEXMK) notes.tex

# Rebuild after each save; open the PDF separately. Stop with Ctrl-C.
watch:
	$(LATEXMK) -pvc -view=none notes.tex

# Remove temporary compilation files; keep the source and finished PDF.
clean:
	$(LATEXMK) -c notes.tex
