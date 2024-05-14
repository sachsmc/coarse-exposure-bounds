# Reproducibility materials for "The impact of coarsening an exposure on partial identifiability in instrumental variable settings" by Erin E Gabriel, Michael C Sachs, and Arvid Sjölander

## Data sources

### Peanut study

Du Toit G, Roberts G, Sayre PH, Bahnson HT, Radulovic S, Santos AF, Brough HA, Phippard D, Basting M, Feeney M, Turcanu V. Randomized trial of peanut consumption in infants at risk for peanut allergy. New England Journal of Medicine. 2015 Feb 26;372(9):803-13. (https://www.nejm.org/doi/full/10.1056/NEJMoa1414850).

The publicly available trial data were downloaded from the Immune Tolerance Network TrialShare website on 2020-06-15 (https://www.itntrialshare.org/, study identifier: ITN032AD). They are not redistributed here, but you can download them after registering at the ITN trialshare website. 

The datasets you will need are ADCONSP2 and ADPOUT1. The analysis dataset is produced from those files in code/01-munge-data.R. 

The analysis is reproduced in the code/03-apply-bounds-peanut.R script. 

### Homocysteine and cardiovascular disease study

Meleady, R., Ueland, P.M., Blom, H., Whitehead, A.S., Refsum, H., Daly, L.E., Vollset, S.E., Donohue, C., Giesendorf, B., Graham, I.M. and Ulvik, A., 2003. Thermolabile methylenetetrahydrofolate reductase, homocysteine, and cardiovascular disease risk: the European Concerted Action Project. The American journal of clinical nutrition, 77(1), pp.63-70 

This paper reports the results of an observational study designed to investigate the effect of homocysteine on cardiovascular disease using the 677CT polymorphism (rs1801133) in the Methylenetetrahydrofolate Reductase gene as an instrument. We extracted the summary data reported in Table 3 of that paper and provide it here in the rawdata subdirectory. 

The analysis is reproduced in the code/03-apply-bounds-mendelian.R.

## Bounds code

Code to compute the bounds are included as R objects saved as rds files in the data subdirectory. They are read in by the scripts automatically. 

The following packages are required: 

causaloptim, 
    here, 
    data.table, 
    ggplot2, 
    ivtools

## License

MIT License, see file LICENSE.md

