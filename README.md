# Digital elevation model quality on digital soil mapping prediction accuracy

## General description

This repository contains the data and R source code used to prepare the paper [__Digital elevation model quality on digital soil mapping prediction accuracy__][paper], published in the journal [Ciência e Agrotecnologia][cagro]. The study was developed by Elias Mendes Costa, Alessandro Samuel-Rosa, and Lúcia Helena Cunha dos Anjos as part of the Master Thesis of Elias Mendes Costa presented before the Post-Graduate Course in Agronomy-Soil Science of the Federal Rural University of Rio de Janeiro on 26 February 2015.

[paper]: http://dx.doi.org/10.1590/1413-70542018426027418
[cagro]: http://www.scielo.br/scielo.php?script=sci_serial&pid=1413-7054&lng=en&nrm=iso

[Elias Mendes Costa](https://www.researchgate.net/profile/Elias_Costa6) was supported by the FAPERJ Foundation
(Process E-26/100.422/2014) and later by the CNPq Foundation (141391/2015-4).
[Alessandro Samuel-Rosa](https://www.researchgate.net/profile/Alessandro_Samuel-Rosa) was supported by the 
CAPES Foundation (Process 88887.116157/2016-00).
[Lúcia Helena Cunha dos Anjos](https://www.researchgate.net/profile/Lucia_Anjos) was supported by the CNPq
Foundation (Process 480515/2013-1).

## Study summary

Digital elevation models (DEM) used in digital soil mapping (DSM) are commonly selected based on measures and indicators (quality criteria) that are thought to reflect how well a given DEM represents the terrain surface. The hypothesis is that the more accurate a DEM, the more accurate will be the DSM predictions. The objective of this study was to assess different criteria to identify the DEM that delivers the most accurate DSM predictions. A set of 10 criteria were used to evaluate the quality of nine DEMs constructed with different data sources, processing routines and three resolutions (5, 20, and 30 m). Multinomial logistic regression models were calibrated using 157 soil observations and terrain attributes derived from each DEM. Soil class predictions were validated using leave-one-out cross-validation. Results showed that, for each resolution, the quality criteria are useful to identify the DEM that more accurately represents the terrain surface. However, for all three resolutions, the most accurate DEM did not produce the most accurate DSM predictions. With the 20-m resolution DEMs, DSM predictions were five percentage points less accurate when using the more accurate DEM. The 5-m resolution was the most accurate DEM overall and resulted in DSM predictions with 44% accuracy; this value was equal to that obtained with two coarser resolution, lower accuracy DEMs. Thus, identifying the truly best DEM for DSM requires assessment of the accuracy of DSM predictions using some form of external validation, because not necessarily the most accurate DEM will produce the best DSM predictions.

## Repository structure

We use the following repository structure:

    dem4dsm
    - code  # R source code
    - data  # point soil data, elevation ground control points, and raster covariate data
    - res   # resulting images

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.
