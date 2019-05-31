# [Disentangling data discrepancies with integrated population models](https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.2714)

### Sarah P. Saunders, Matthew T. Farr, Alexander D. Wright, Christie A. Bahlai, José W. Ribeiro Jr., Sam Rossman, Allison L. Sussman, Todd W. Arnold, & Elise F. Zipkin

### Ecology (*In Press*)

### Code/Data DOI: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2532006.svg)](https://doi.org/10.5281/zenodo.2532006)

### Please contact the first author for questions about the code or data: Sarah P. Saunders (sarahpsaunders@gmail.com)
__________________________________________________________________________________________________________________________________________

## Abstract:
A common challenge for studying wildlife populations occurs when different observation methods (e.g. mark-recapture, census, telemetry, occupancy) provide inconsistent or incomplete inference on the trend, dynamics, or viability of a population. A potential solution to the challenge of conflicting or piecemeal data (i.e. dataset discrepancies and deficiencies, respectively) relies on the integration of multiple data types into a unified modeling framework, such as integrated population models (IPMs). IPMs are a powerful approach for species that inhabit spatially and seasonally complex environments. In cases where data come from heterogeneous collection schemes, individual data sources are unlikely to contain all the information necessary to model large-scale population dynamics, and may fail to capture spatial and temporal variability across a species’ range or annual cycle. We review how data discrepancies can be modeled with IPMs and provide guidance on combining multiple data types when there are data deficiencies. We illustrate this issue with analysis of a migratory species, the American woodcock (Scolopax minor), in which individual monitoring programs suggest differing population trends. To address this discrepancy, we synthesized several long-term datasets (1963 – 2015) within an IPM to estimate continental-scale population trends, and link dynamic drivers across the full annual cycle and complete extent of the woodcock’s geographic range in eastern North America. Our analyses reveal the limiting portions of the life cycle by identifying time periods (i.e. seasons) and locations (i.e. population units) where survival rates are lowest and most variable, as well as which demographic parameters constitute the main drivers of population change. We provide recommendations for exploiting the capabilities of IPMs to account for incongruent datasets and discuss how strategies (e.g. data weighting, data thinning) from other disciplines can be incorporated into ecological analyses that combine multiple, disparate data types.  

## Code 
1. [Data_manip_AMWO_marray_Git.R](https://github.com/zipkinlab/timberdoodle/blob/master/Data_manip_AMWO_marray_Git.R): Data managment and manipulation code of Banding data to create m-arrays for band-recovery model component of the IPM. 
2. [AMWO_harvest_Git.R](https://github.com/zipkinlab/timberdoodle/blob/master/AMWO_harvest_Git.R): Data manipulation of harvest data (Duck Stamp Survey, Harvest Information Program, and Parts-collection Surveys) and JAGS model file and R code for estimation of annual woodcock harvest throughought temporal and spatial range of the study.
3. [IPM_AMWO_modelfile_Git.R](https://github.com/zipkinlab/timberdoodle/blob/master/IPM_AMWO_modelfile_Git.R): JAGS model file and R code for bayesian implementation of the Integrated Population model to estimate population dynamics and trends of the Eastern and Central US woodcock populations from 1963-2015.

## Data
This project utilizes 3 datasets:  

1) Banding-Recovery Data  
Description: Bird banding is a universal technique for studying the movement, survival and behavior of birds. Recovery data is also useful for estimating survival rates of both young and adult birds.  
Contact: USGS Bird Banding Lab  
Website: https://www.pwrc.usgs.gov/bbl/  

2) Hunter Harvest Data  
Description: The U.S. Harvest Information Program (HIP) measures hunter effort and number of woodcocks harvested by hunters annually (1963-Present).  
Contact: USFWS Migratory Bird Data Center  
Website: https://migbirdapps.fws.gov/  

3) Singing Ground Survey  
Description: Survey of numerous routes throughout the range of woodcocks to survey singing males during the spring.   	
Contact: USFWS Migratory Bird Data Center  
Website: https://migbirdapps.fws.gov/  

## Raw data files

1) [AMWO bandings](https://github.com/zipkinlab/Saunders_etal_2019_Ecol/blob/master/Raw_Data/AMWO%20bandings.csv): The USGS Bird Banding Laboratory has compiled banding data for woodcocks over the complete timeframe of our study period (1963 – 2015).

2) [AMWO recoveries](https://github.com/zipkinlab/Saunders_etal_2019_Ecol/blob/master/Raw_Data/AMWO%20recoveries.csv): The USGS Bird Banding Laboratory has compiled recovery data for woodcocks over the complete timeframe of our study period (1963 – 2015).

3) [Harvest model data](https://github.com/zipkinlab/Saunders_etal_2019_Ecol/blob/master/Raw_Data/Harvest_model_data.Rda): Model-based estimates from our harvest model using Duck Stamp Survey (1963 - 2001) and Harvest Information Program data (1999 - 2015) on American woodcock.

4) [SGS indices](https://github.com/zipkinlab/Saunders_etal_2019_Ecol/blob/master/Raw_Data/SGS-indices.csv): Singing-ground Survey indices of American woodcock during 1968 - 2015 from US Fish and Wildlife Service,

## Data files produced from models and referred to in the code:

1) AMWO marray: a multidimensional band-recovery array (referred to as an ‘m-array’) consisting of 12 subcomponent arrays

2) relAMWO: a separate array of total birds banded during each year, banding period, stage-sex class, and population unit, which was provided as data to the model

3) AMWO harvest pi Oct 2017: the saved model output from the full IPM model run, used to create figures and results visualizations

4) [AMWO pop sizes](https://github.com/zipkinlab/Saunders_etal_2019_Ecol/blob/master/AMWO_pop_sizes.csv): Annual population size estimates from our IPM for both Central and Eastern management units during 1963 - 2015. Note that estimates are in the millions (the scale of estimates shown in Figure 2).

