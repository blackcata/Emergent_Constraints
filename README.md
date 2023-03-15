# Emergent Constraints
Basic python code for calculating the emergent constraint method, which is based on the Cox et al., 2013. 
The emergent constraint method is reducing the uncertainty or diversity of future projections in individual models by current states (Hall et al., 2019).
Due to diverse equations and parameterizations, the errors are accumulated through time integration. 
Therefore, the Earth system models (ESMs) often projects diversity even in a single model, and larege diversity in inter-model dimension.
Then, we cannot trust the future projections in ESMs. We need to reduce the uncertainty of future projections to make reliable ranges.
The emergent constrain method can reduce the uncertainty of future projections by relationship between current state and future changes.
For example, our paper represents the relationship between current nitrate climatology and future Arctic chlorophyll and productivity changes.

This code's Arctic chlorophyll estimattions are published at the Earth's Future (AGU)

### Title : Emergent Constraint for Future Decline in Arctic Phytoplankton Concentration

Listed co-author(s) : Kyung Min Noh, Hyung-Gyu Lim, Eun-Jin Yang, and Jong-Seong Kug

Corresponding Author : Hyung-Gyu Lim, and Jong-Seong Kug
  
### 1. Provided Data 
CHL, NO3, PP projections from the CMIP5 and CMIP6 were used for emergent constrained method.
Additionally, the climatology of NO3 was used to constrain the estimations in CHL projections.
You can download these reanalysis and CMIP data from below URL.
  - CMIP     : https://esgf-node.llnl.gov/
  - WOA18    : https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/bin/woa18oxnu.pl?parameter=n 
  - GLODAPv2 : https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0162565/


### 2. Code composition
  - emergent_constraint.py [basic python functions for emergent constrained method]
  - Emergent_Constraints.ipynb [Jupyter codes for drawing plots]
  
Before start the codes, you need to download the dataset to DATA folder, which is located same path on python and jupyter code.
 
All output pdf files will be saved on the RESULT directory.
