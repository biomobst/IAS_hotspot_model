# IAS Hotspot Model

## Summary 

<img src=images/frontpage3.jpg width=500 align=right>

The scripts in this repo are designed to use species distribution models (SDM) to identify potential high-risk areas for the introduction, establishment and spread of invasive species in Swedish water. The SDM workflow is based on a modeling workflows developed by University of Gothenburg (Leidenberger et al. 2015; Laugen et al. 2015; Stelzer et al. 2013; Karlsson et al. 2019) and was further adapted to model a large number of known invasive species that have not yet established themselves in Sweden, so-called "door knocking species". The models' results not only identify potential areas of distribution for each individual species, but can also be used to map regions where suitable habitats for a large number of invasive species overlap. These superimposed maps are then further integrated with external data layers on introduction vectors such as ship traffic density and oceanographic currents as well as with data on invasive impact for individual species. The integration of oceanographic data is however is not yet fully implemented. The regions identified in this workflow as areas with a high overall invasion risk (i.e. high risk of both introduction and establishment) can be considered invasive hotspots.

The workflow for running the analysis is written in R and was built under the commission of the Swedish Agency for Marine and Water Management (SwAM). It is formally part of the report **“Provtagningsdesign med omdrev för övervakning av främmande arter enligt eRAS”** by Marine Monitoring AB, Senalytics AB, and the Swedish National Veterinary Institute. The repo also contains all input and output data from the original experiments. Here is the link to the report: https://www.havochvatten.se/data-kartor-och-rapporter/rapporter-och-andra-publikationer/publikationer/2020-10-09-provtagningsdesign-for-overvakning-av-frammande-arter.html

We also customised of the original IAS hotspot model to predict suitable habitat in the northern Baltic Sea for invasive species originating from freshwater systems, i.e. lakes and rivers. This model is documented in the folder "Freshwater model for Bothnian sea" above.

## References

Karlsson R, Obst M, Berggren M (2019) Analysis of potential distribution and impacts for two species of alien crabs in Northern Europé Biological Invasions. https://link.springer.com/article/10.1007/s10530-019-02044-3

Laugen AT, Hollander J, Obst M, Strand A (2015) The Pacific Oyster (Crassostrea gigas) invasion in Scandinavian coastal waters in a changing climate: impact on local ecosystem services. In Biological Invasions in Aquatic and Terrestrial Systems: Biogeography, Ecological Impacts, Predictions, and Management. De Gruyter, Warsaw. Pp. 232-257. https://www.researchgate.net/publication/290061791_The_Pacific_Oyster_Crassostrea_gigas_invasion_in_Scandinavian_coastal_waters_impact_on_local_ecosystem_services

Leidenberger S, Obst M, Kulawik R, Stelzer K, Heyer K, Hardisty A, Bourlat SJ (2015) Evaluating the potential of ecological niche modelling as a component in non-indigenous species risk assessments. Marine Pollution Bulletin. 97: 470-487. https://www.sciencedirect.com/science/article/pii/S0025326X15002350

Stelzer K, Heyer K, Bourlat S, Obst M (2013) Application of Niche Modeling and Earth Observation for the risk assessment and monitoring of invasive species in the Baltic Sea. Report MarCoast II - Marine and Coastal Environmental Information Services, Ballast Water Option, pp 57. https://zenodo.org/record/886349#.XvCfe1BS_fZ

Strand M, Aronsson M, Svensson M (2018) Klassificering av främmande arters effekter på biologisk mångfald i Sverige – ArtDatabankens risklista. ArtDatabanken Rapporterar 21. ArtDatabanken SLU, Uppsala. Hämtad 2020-05-04 från file://hav.havochvatten.se/hav/root/users/erllet/Documents/REF/Strand%20m%20fl_Riskklistan%20ADb%202018.pdf
