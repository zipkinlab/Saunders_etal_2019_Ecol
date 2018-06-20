## [Disentangling data discrepancies and deficiencies with integrated population models](https://msu.edu/user/ezipkin/)  
#### *Sarah P. Saunders, Matthew T. Farr, Alexander D. Wright, Christie A Bahlai, Jose W. Ribeiro Jr., Sam Rossman, Allison L. Sussman, Todd W. Arnold, & Elise F. Zipkin*  
*Permanent Identifier (e.g. DOI):* ...  
**Abstract:** A common challenge for studying wildlife populations occurs when different observation methods (e.g. mark-recapture, census, telemetry, occupancy) provide inconsistent or incomplete inference on the trend, dynamics, or viability of a population. A potential solution to the challenge of conflicting or piecemeal data (i.e. dataset discrepancies and deficiencies, respectively) relies on the integration of multiple data types into a unified modeling framework, such as integrated population models (IPMs). IPMs are a powerful approach for species that inhabit spatially and seasonally complex environments. In cases where data come from heterogeneous collection schemes, individual data sources are unlikely to contain all the information necessary to model large-scale population dynamics, and may fail to capture spatial and temporal variability across a species’ range or annual cycle. We review how data discrepancies can be modeled with IPMs and provide guidance on combining multiple data types when there are data deficiencies. We illustrate this issue with analysis of a migratory species, the American woodcock (Scolopax minor), in which individual monitoring programs suggest differing population trends. To address this discrepancy, we synthesized several long-term datasets (1963 – 2015) within an IPM to estimate continental-scale population trends, and link dynamic drivers across the full annual cycle and complete extent of the woodcock’s geographic range in eastern North America. Our analyses reveal the limiting portions of the life cycle by identifying time periods (i.e. seasons) and locations (i.e. population units) where survival rates are lowest and most variable, as well as which demographic parameters constitute the main drivers of population change. We provide recommendations for exploiting the capabilities of IPMs to account for incongruent datasets and discuss how strategies (e.g. data weighting, data thinning) from other disciplines can be incorporated into ecological analyses that combine multiple, disparate data types.


### Code 
File Name | File Location | File Description
--- | --- | ---
Data_manip_AMWO_marray_Git.R 	| _Final submission files June 2018_ folder  |  Data managment and manipulation code of Banding data to create m-arrays for band-recovery model component of the IPM. 
AMWO_harvest_Git.R 		| _Final submission files June 2018_ folder  |  Data manipulation of harvest data (Duck Stamp Survey, Harvest Information Program, and Parts-collection Surveys) and JAGS model file and R code for estimation of annual woodcock harvest throughought temporal and spatial range of the study.
IPM_AMWO_modelfile_Git.R 	| _Final submission files June 2018_ folder  | JAGS model file for bayesian implementation of the Integrated Population model to estimate population dynamics and trends of the Eastern and Central US woodcock populations from 1963-2015.

### Data
This project utilizes 3 datasets:  

1) Banding-Recovery Data  
&nbsp;&nbsp;Description: Bird banding is a universal technique for studying the movement, survival and behavior of birds. Recovery data is also useful for estimating survival rates of both young and adult birds.  
&nbsp;&nbsp;Contact: USGS Bird Banding Lab  
&nbsp;&nbsp;Website: https://www.pwrc.usgs.gov/bbl/

2) Hunter Harvest Data  
&nbsp;&nbsp;Description: The U.S. Harvest Information Program (HIP) measures hunter effort and number of woodcocks harvested by hunters annually (1963-Present).  
&nbsp;&nbsp;Contact: USFWS Migratory Bird Data Center  
&nbsp;&nbsp;Website: https://migbirdapps.fws.gov/

3) Singing Ground Survey  
&nbsp;&nbsp;Description: Survey of numerous routes throughout the range of woodcocks to survey singing males during the spring.     	
&nbsp;&nbsp;Contact: USFWS Migratory Bird Data Center  
&nbsp;&nbsp;Website: https://migbirdapps.fws.gov/

### Contact
Funded by University of Minnesota, Migratory Shore and Upland Game Bird Grant (USFWS Webless Migratory Game Bird Research and Management Program).  
&nbsp;&nbsp;&nbsp;&nbsp;PI: Dr. Todd Arnold  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Department of Fisheries, Wildlife, and Conservation Biology  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;University of Minnesota  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;St. Paul, MN 55108  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;612-624-2220  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;arnol065@umn.edu  
