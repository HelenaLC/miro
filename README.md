# miro

[DEMO](https://helenalc.github.io/miro/articles/miro.html)

`miro` provides a simple and effective approach for EDA 
based on color representations of molecular variation via: 

1. Rescaling transcription-derived principal components (PCs) 
so they can be interpreted as colors in a perceptually meaningful way.
2. Visualizing these colors in Euclidean space, yielding a tissue image.
3. Rotations of the underlying PC space help identify which 
features contribute most strongly to observed spatial patterns.

In the **CIELAB color space** (or $L^*a^*b^*$), distances between colors 
approximately reflect differences perceived by the human eye. 
Each color is represented by three components: $L$ for lightness, as well as
$a$ and $b$ for the red–green and yellow–blue color axes, respectively.

![**Visium data on colorectal cancer.**
**(a)** Hematoxylin and eosin (H&E) staining, full section (top) and ST
region (bottom); epi. = intestinal epithelium, tum. = colorectal cancer. 
**(b)** LAB colors are expressed as $a\in[-100,100]$ (green–red), 
$b\in[-100,100]$ (blue–yellow), and $L\in[0,100]$ (dark–light).
**(c)** Scatter plot of principal components (PCs) 1–3 with points (= cells) colored by LAB.
**(d)** Spatial plots with spots colored by the LAB interpretation of PCs 1-3 as ABL, and PCs 1-3.
**(e)** Biplot from the LAB interpretation of feature loadings for PCs 1-2 as AB,
and L fixed to 75. Included are selected canonical markers in this context.](vignettes/schematic.png)
