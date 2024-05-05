# IAS Hotspot Model

## Summary 

<img src=images/frontpage3.jpg width=500 align=right>

The scripts in this repo are designed to use species distribution models (SDM) to identify potential high-risk areas for the introduction, establishment and spread of invasive species in European coastal water. The SDM workflow is based on a modeling workflow developed by University of Gothenburg (Leidenberger et al. 2015; Laugen et al. 2015; Stelzer et al. 2013; Karlsson et al. 2019) and was further adapted to model a large number of known invasive species. The models' results not only identify potential areas of distribution for individual species, but can also be used to map regions where suitable habitats for a large number of invasive species overlap. These superimposed maps can then further be integrated with external data layers on introduction vectors such as ship traffic density and oceanographic currents as well as with data on invasive impact for individual species. The integration of oceanographic data is however is not yet fully implemented. The regions identified in this workflow as areas with a high overall invasion risk (i.e. high risk of both introduction and establishment) can be considered invasive hotspots.

The workflow for running the analysis is written in R and was built under the commission of the Swedish Agency for Marine and Water Management (SwAM) and was futher supported by the European project MARCO-BOLO (MARine COastal BiOdiversity Long-term Observations). It is formally part of the reports by Bergkvist et al (2020) and Obst & Andersson (2023). The repo contains all input and output data from the modelling experiments and is continuously developed further.

The original IAS hotspot model was expanded to predict suitable habitat in the northern Baltic Sea for invasive species originating from freshwater systems, i.e. lakes and rivers. This model is documented in the folder "Freshwater model for Bothnian sea".

## References

Obst M, Andersson G (2023) Hotspot modell för invasiva arter i Bottniska viken - Övervakning i marin miljö. SeAnalytics rapport 2023-03. På uppdrag av Havs- och vattenmyndigheten.

Bergkvist J, Magnusson M, Obst M, Sundberg P, Andersson G (2020) Provtagningsdesign för övervakning av främmande arter. Övervakning i marin miljö. Havs- och vattenmyndighetens rapport 2020:22. ISBN 978-91-88727-86-2

Karlsson R, Obst M, Berggren M (2019) Analysis of potential distribution and impacts for two species of alien crabs in Northern Europé Biological Invasions. https://link.springer.com/article/10.1007/s10530-019-02044-3

Laugen AT, Hollander J, Obst M, Strand A (2015) The Pacific Oyster (Crassostrea gigas) invasion in Scandinavian coastal waters in a changing climate: impact on local ecosystem services. In Biological Invasions in Aquatic and Terrestrial Systems: Biogeography, Ecological Impacts, Predictions, and Management. De Gruyter, Warsaw. Pp. 232-257. https://www.researchgate.net/publication/290061791_The_Pacific_Oyster_Crassostrea_gigas_invasion_in_Scandinavian_coastal_waters_impact_on_local_ecosystem_services

Leidenberger S, Obst M, Kulawik R, Stelzer K, Heyer K, Hardisty A, Bourlat SJ (2015) Evaluating the potential of ecological niche modelling as a component in non-indigenous species risk assessments. Marine Pollution Bulletin. 97: 470-487. https://www.sciencedirect.com/science/article/pii/S0025326X15002350

Stelzer K, Heyer K, Bourlat S, Obst M (2013) Application of Niche Modeling and Earth Observation for the risk assessment and monitoring of invasive species in the Baltic Sea. Report MarCoast II - Marine and Coastal Environmental Information Services, Ballast Water Option, pp 57. https://zenodo.org/record/886349#.XvCfe1BS_fZ

Strand M, Aronsson M, Svensson M (2018) Klassificering av främmande arters effekter på biologisk mångfald i Sverige – ArtDatabankens risklista. ArtDatabanken Rapporterar 21. ArtDatabanken SLU, Uppsala. Hämtad 2020-05-04 från file://hav.havochvatten.se/hav/root/users/erllet/Documents/REF/Strand%20m%20fl_Riskklistan%20ADb%202018.pdf
