# Seasonal constraints on inferred planetary heat content
McKinnon and Huybers (2016), GRL

This repo contains the codes (and some limited data) used for the calculations and figures in McKinnon and Huybers (2016).

The codes are in ncl, matlab, and python. There is only one matlab code which can be easily skipped if you do not have access to matlab by simply loading the argo.missing.mask.csv (in the data folder) rather than running the matlab code.

The first step is to run the ncl scripts that are wrapped using calcClimos.ncl. This script processes most of the data and translates the physical quantities into heat units. The original data sources are all referenced in the relevant ncl scripts.

The remainder of the analysis, including the creation of the figures for the paper, is performed via an ipython notebook, compareSeasonalCycles. 

As always, I may have forgotten a subfunction or two, so please be in touch if something is missing. 

Questions/comments/concerns/ideas: mckinnon@ucar.edu
