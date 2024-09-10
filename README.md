# Short-term effect of air pollution on OHCA in Lombardia – a case-crossover study
It is time-stratified case-crossover study to analyse the association between short-term exposure to air pollution and Out-of-Hospital-Cardiac-Arrests (OHCA) between 2016 and 2019 in the region of Lombardia, Italy. 
The main model consists of two-stages of multi-scaler risk assessment:
## Stage 1:
Estimating OHCA risk for 96 districts (of approximately 100k population) using a distributed lag non-linear model (DLNM) with individual air pollutant as exposure factor, and temperature and relative humidity as confounders. Long-term trends and seasonality are controlled for using a natural cubic spline of time. 
## Stage 2:
Estimates from ditricts are pooled in to calculate the effect for the entire region using random-effect meta-analysis.
< br / > 
Outputs are produced as:
* exposure-response function illustrating changing risk over increasing exposure to the pollutant,
* lag-response function illustrating changing risk per 10 unit increase in the pollutant exposure between lag 0 (event day) and lag 7.
   
## Data
### Air pollutants 
* PM2.5, PM10, NO2, O3, SO2 and CO (in µg/m3)
* Spatial resolution - 0.1° x 0.1° (approximately 10 km x 10 km)
* Temporal resolution - 24 hours
### OHCA 
Geolocalised points of telephone calls to the regional emergency medical service made for cardiac arrests. 

## Files information
1. Data : contains aggregated district-wise data of daily air pollution concentration, temperature, relative humidity, OHCA events, date and case-crossover spectrum 
2. Graphs and maps : outputs as dsitrict-wise relative risk maps, regional exposure-response curves and regional lag-response graphs
3. Results: excel consisiting relative risk (ci:95%) (for entire region) for overall model, stratified by urban-rural, modified by age and sex and bi-pollutant. 
4. Parametres_selection.Rmd : Selections of model parameters such as degrees of freedom, lag and spline function of temperature, relative humidity and pollutants.
5. Model: consists R code for 
> * Main model
> * Model modification by age-sex
> * Model stratification by urban-rural
> * Bi-pollutant model 
