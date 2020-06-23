# IAS Hotspot Model

## Summary 

<img src=images/frontpage3.jpg width=400 align=right>

This repo contains a workflow for running species distribution models for marine species in R and use it to identify potential distribution hotspots areas for alien / invasive species in Swedish waters. The repo also contains all input and output data from the original experiments.

The IAS Hotspot Model is based on Species Distribution Models (SDM) embedded in a workflow that allows predicting suitable habitat for a large number of invasive species. These individual predictions are summarized and can subsequently be combined with external data layers on introduction vectors such as ship traffic density and oceanographic currents. The latter option however is not yet fully implemented. 

SDMs are based on a machine learning algorithm (Random Forest) and analyze the correlation between environmental variables and species occurrences. Predictions are carried out with the same environmental variables as during the model building, a so-called "native projections". A more detailed account of the modeling theory can be found in the below papers.

The workflow was built under the commission of the Swedish Agency for Marine and Water Management (SwAM)” and is formally part of the report **“Provtagningsdesign med omdrev för övervakning av främmande arter enligt eRAS”** by Marine Monitoring AB, Senalytics AB, and the Swedish National Veterinary Institute.

## References

Karlsson R, Obst M, Berggren M (2019) Analysis of potential distribution and impacts for two species of alien crabs in Northern Europé Biological Invasions. https://link.springer.com/article/10.1007/s10530-019-02044-3

Leidenberger S, Obst M, Kulawik R, Stelzer K, Heyer K, Hardisty A, Bourlat SJ (2015) Evaluating the potential of ecological niche modelling as a component in non-indigenous species risk assessments. Marine Pollution Bulletin. 97: 470-487. https://www.sciencedirect.com/science/article/pii/S0025326X15002350

Stelzer K, Heyer K, Bourlat S, Obst M (2013) Application of Niche Modeling and Earth Observation for the risk assessment and monitoring of invasive species in the Baltic Sea. Report MarCoast II - Marine and Coastal Environmental Information Services, Ballast Water Option, pp 57. https://zenodo.org/record/886349#.XvCfe1BS_fZ

Laugen AT, Hollander J, Obst M, Strand A (2015) The Pacific Oyster (Crassostrea gigas) invasion in Scandinavian coastal waters in a changing climate: impact on local ecosystem services. In Biological Invasions in Aquatic and Terrestrial Systems: Biogeography, Ecological Impacts, Predictions, and Management. De Gruyter, Warsaw. Pp. 232-257. https://www.researchgate.net/publication/290061791_The_Pacific_Oyster_Crassostrea_gigas_invasion_in_Scandinavian_coastal_waters_impact_on_local_ecosystem_services
