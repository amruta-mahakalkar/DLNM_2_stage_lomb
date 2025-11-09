# Short-term effect of air pollution on OHCA in Lombardia – a case-crossover study

It is time-stratified case-crossover study to analyse the association between short-term exposure to air pollution and Out-of-Hospital-Cardiac-Arrests (OHCA) between 2016 and 2019 in the region of Lombardia, Italy published thorugh an article at https://doi.org/10.1002/gch2.202500241
The main model consists of two-stages of multi-scaler risk assessment:

### Stage 1

Estimating OHCA risk for 96 districts (of approximately 100k population) using a [distributed lag non-linear model](https://github.com/gasparrini/dlnm) (DLNM) with individual air pollutant as exposure factor, and temperature and relative humidity as confounders. Long-term trends and seasonality are controlled for using a stratum (yyyy:month:dow) as per the case-crossover study design.

### Stage 2

Estimates from districts are pooled in to calculate the effect for the entire region using random-effect meta-analysis.

Outputs are produced as:

* Exposure-response function illustrating changing risk over increasing exposure to the pollutant,
* Lag-response function illustrating changing risk per 10 unit increase in the PM2.5, PM10, NO2 and O3, 1 unit increase in SO2 and 100 unit increase in CO exposure between lag 0 (event day) and lag 7.

## Data

### Exposure

* PM2.5, PM10, NO2, O3, SO2 and CO (in µg/m3)
* Spatial resolution - 0.1° x 0.1° (approximately 10 km x 10 km)
* Temporal resolution - 24 hours

### OHCA

Geolocalised points of telephone calls to the regional emergency medical service made for cardiac arrests.

## Files information

1. Graphs and maps : outputs as district-wise relative risk maps, regional exposure-response curves and regional lag-response graphs
2. Models and other: contains R codes for

   * Data preparation
   * Sensitivity analysis for parameter selection
   * Main model
   * Model stratification by age-sex
   * Model stratification by urban-rural
   * Model stratification by seasons
   * Bi-pollutant model
   * Z-test to check effect modification by age, sex, location and seasons.
   * Codes for plotting exposure-lag-response relationships
