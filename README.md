# PIJA_Management
Data and code supporting Steel et al. "Range-wide effects of vegetation management on a rapidly declining woodland bird"

# Manuscript Abstract: 
Birds face population declines globally from a combination of threats including habitat loss and degradation. In woodlands and forests of the western U.S., vegetation management is expanding to mitigate wildfire risk and meet other objectives as historic droughts increase ecological stressors and conservation challenges. These intensifying threats create an urgent need to evaluate the impacts of multiple vegetation management techniques and how those effects vary across time and habitat context. We demonstrate a pathway to address these information needs using a case study for the declining pinyon jay (*Gymnorhinus cyanocephalus*), a species petitioned for listing under the U.S. Endangered Species Act. We combined 52,702 participatory science bird checklists (2011–2022) with spatially explicit treatment data to quantify effects of prescribed fire, mastication, thinning, and timber harvest on pinyon jay counts across the species' U.S. range. We then evaluated how treatment effects varied by season, time since treatment, and baseline density of key tree species. During the breeding season, jays generally avoided treated landscapes, with prescribed fire and thinning showing clear negative effects. Impacts were more varied during the non-breeding season, with mastication and recent prescribed fires showing positive associations with jay counts in sparse woodlands, while thinning and harvest showed consistent negative effects. These findings provide timely information at broad scales for the management of a species of conservation concern and its habitat. Further, our treatment-comparative, context-dependent approach offers a transferable method to evaluate biodiversity tradeoffs when multiple vegetation objectives are pursued across heterogeneous landscapes. 

# Code notes:
This repository contains data prepped for running the two main zero-inflated Poisson models used in the manuscript as well as code to run and interpret model outputs, including manuscript figures.  

Rerunning the models should give you *nearly* identical results. The stochastic MCMC process may generate very minor differences, but we have also replaced the original CCI observer scores with the mean value of the metric. This is required to protect the privacy of eBird participants. Calculations for the checklist calibration index are described in Kelling et al. 2015 and Johnston et al. 2018. 

Kelling, S., A. Johnston, W. M. Hochachka, M. Iliff, D. Fink, J. Gerbracht, C. Lagoze, F. A. L. Sorte, T. Moore, A. Wiggins, W.-K. Wong, et al. (2015). Can Observation Skills of Citizen Scientists Be Estimated Using Species Accumulation Curves? PLOS ONE 10:e0139600.

Johnston, A., D. Fink, W. M. Hochachka, and S. Kelling (2018). Estimates of observer expertise improve species distributions from citizen science data. Methods in Ecology and Evolution 9:88–97.

# How to cite:
If the data or code associated with this repository are used, please cite the associated manuscript:

Zachary L. Steel, Hailey M. Boone, Leah McTigue, Valerie Stein Foster, Merijn van den Bosch, Andrew N. Stillman, Caroline D. Cappello, Mark Ditmer, and Emily J. Francis. "Range-wide effects of vegetation management on a rapidly declining woodland bird". (Currently in review)
