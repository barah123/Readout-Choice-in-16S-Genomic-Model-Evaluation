# PLOS Computational Biology submission figures

Seven figures, named as the journal expects. Regenerate with
`python ../../scripts/make_figures.py`.

| file | content | width (in) | height (in) |
|---|---|---|---|
| Fig1.tiff | Paired differences, locked test and replications | 7.20 | 4.60 |
| Fig2.tiff | Readout, model and layer effects, metrics kept separate | 7.20 | 6.20 |
| Fig3.tiff | All methods, bounded by permuted labels and the ceiling | 4.60 | 3.40 |
| Fig4.tiff | Task difficulty and the hard-panel extension | 7.00 | 3.20 |
| Fig5.tiff | Sensitivity to sequence degradation | 5.40 | 3.30 |
| Fig6.tiff | Region-selection regret | 7.20 | 3.30 |
| Fig7.tiff | Achieved accuracy against the ambiguity ceiling | 4.80 | 4.20 |

## Compliance

- **Format** TIFF, LZW compression
- **Resolution** 300 dpi (PLOS accepts 300 to 600)
- **Colour** RGB, 8-bit, alpha flattened onto white
- **Width** all within the 2.63 to 7.5 inch range
- **Height** all well under the 8.75 inch maximum
- **Size** 1.3 MB total, every file far below the 10 MB limit
- **Text** 7 to 10 pt at final size
- **Palette** Okabe-Ito, colourblind-safe; no information is carried by colour alone in
  Figs 1, 4, 5 and 6, where position or label also distinguishes the series

## Two notes for production

1. Figure 2 was reflowed from a single row of four panels to a two by two grid. The single
   row measured 9.6 inches, above the PLOS maximum, and simply scaling it down would have
   pushed labels below the 8 pt minimum.
2. Figures are saved via an RGBA render that is then flattened onto white. Writing TIFF
   directly from matplotlib leaves an alpha channel, which PLOS production does not expect.

Working copies in PNG and vector PDF are in the parent directory.
