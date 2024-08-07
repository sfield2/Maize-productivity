# Maize-productivity
All scripts relating to Field et al. 2024 JCAA Submission.
This directory contains code for:
1. Running the maize paleoproductivity model.
2. Calculating maize surpluses or deficits depending on community needs.
3. Running summary analyses and building figures.

# Data
Portions of the data used in this research is provided infolder entitled "DATA". Sensitive site location data is not included.
Spatial data  (catchment.shp, drainage_area_high.shp, drainage_area_moderate.shp, drinage_area_low.shp, and mesa_top_area.shp) do not represent actual locations of farming areas and catchments used in the study, due to site sensitivity data. Instead, these data are "invented" to show how the model functioned and were are located in publicly accessible areas of the Far View community in Mesa Verde National Park. Climate data used in these steps were downloaded at https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/northamerica/usa/bocinsky2016/ and are derived from Paleocar (https://github.com/bocinsky/paleocar). All tabular data, including complete results from the model, can be found in the DATA folder and associated subfolders. These data contain no sensitive information. 

# Analysis 
Script used for analysis is an R script, entitled "ANALYSIS.R", that contains all necessary commands for running the maize paleoproductivity model, calculating surpluses, and running summary analyses for figures. 

See [Glowacki and Field (2023)](https://github.com/sfield2/NSJ-MV-CeramicSeriation) for methods related to ceramic-based occupation assignments.


If you have questions, email Sean Field (Sean.Field@uwyo.edu).
