# Disentangling data discrepancies with integrated population models

### Sarah P. Saunders, Matthew T. Farr, Alexander D. Wright, Christie A. Bahlai, José W. Ribeiro Jr., Sam Rossman, Allison L. Sussman, Todd W. Arnold, & Elise F. Zipkin

### Ecology

### Please contact the first author for questions about the code or data: Sarah P. Saunders (sarahpsaunders@gmail.com)
__________________________________________________________________________________________________________________________________________

## Abstract:
A common challenge for studying wildlife populations occurs when different observation methods (e.g. mark-recapture, census, telemetry, occupancy) provide inconsistent or incomplete inference on the trend, dynamics, or viability of a population. A potential solution to the challenge of conflicting or piecemeal data (i.e. dataset discrepancies and deficiencies, respectively) relies on the integration of multiple data types into a unified modeling framework, such as integrated population models (IPMs). IPMs are a powerful approach for species that inhabit spatially and seasonally complex environments. In cases where data come from heterogeneous collection schemes, individual data sources are unlikely to contain all the information necessary to model large-scale population dynamics, and may fail to capture spatial and temporal variability across a species’ range or annual cycle. We review how data discrepancies can be modeled with IPMs and provide guidance on combining multiple data types when there are data deficiencies. We illustrate this issue with analysis of a migratory species, the American woodcock (Scolopax minor), in which individual monitoring programs suggest differing population trends. To address this discrepancy, we synthesized several long-term datasets (1963 – 2015) within an IPM to estimate continental-scale population trends, and link dynamic drivers across the full annual cycle and complete extent of the woodcock’s geographic range in eastern North America. Our analyses reveal the limiting portions of the life cycle by identifying time periods (i.e. seasons) and locations (i.e. population units) where survival rates are lowest and most variable, as well as which demographic parameters constitute the main drivers of population change. We provide recommendations for exploiting the capabilities of IPMs to account for incongruent datasets and discuss how strategies (e.g. data weighting, data thinning) from other disciplines can be incorporated into ecological analyses that combine multiple, disparate data types.  

## Code 
1. [Data_manip_AMWO_marray_Git.R](): Data managment and manipulation code of Banding data to create m-arrays for band-recovery model component of the IPM. 
2. [AMWO_harvest_Git.R](): Data manipulation of harvest data (Duck Stamp Survey, Harvest Information Program, and Parts-collection Surveys) and JAGS model file and R code for estimation of annual woodcock harvest throughought temporal and spatial range of the study.
3. [IPM_AMWO_modelfile_Git.R](): JAGS model file and R code for bayesian implementation of the Integrated Population model to estimate population dynamics and trends of the Eastern and Central US woodcock populations from 1963-2015.

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
