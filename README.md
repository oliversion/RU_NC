# RU_NC
Natural computing

For this project, our aim is to model an epidemic disease spread simulator that takes into account external factors, such as density of the population and the vaccination rate. Based on the current anti- vaccination movement [2], it is interesting to investigate what are the results of reduction of vaccination rates in highly populated areas. We plan to build a density-map of the Netherlands for this simulation.

Traditional mathematical epidemic simulation models based on differential equations have drawbacks covered in [4]. In order to deal with this problem, we will use 2D cellular automata to simulate healthy/infected population. We will use [4], where the basic CA model for disease spread is explained, as a base for our implementation.

To be able to simulate larger populations, we will use CA with population density explained in [1]. In this paper, the authors used map-based cellular automata where each cell represents a population cluster with a size determined by the density living in that area. This allows to adjust spreading rate to actual amount of citizen in a cluster. After that we plan to combine the proposed model with information about vaccination (vaccination map or vaccination rate depending on a region). The approach to use the simulating disease spread with vaccination is described in [3].

**References**

[1] A. Holko, M. M ̧edrek, Z. Pastuszak, and K. Phusavat. Epidemiological modeling with a population density map-based cellular automata simulation system. Expert Systems with Applications, 48:1 – 8, 2016.

[2] M J Knol, A T Urbanus, E M Swart, L Mollema, W L Ruijs, R S van Binnendijk, M J te Wierik, H E de Melker, A Timen, and S J Hahn ́e. Large ongoing measles outbreak in a religious community in the netherlands since may 2013. Eurosurveillance, 18(36), 2013.

[3] G.Ch Sirakoulis, I Karafyllidis, and A Thanailakis. A cellular automaton model for the effects of population movement and vaccination on epidemic propagation. Ecological Modelling, 133(3):209 – 223, 2000.

[4] S. Hoya White, A. Mart ́ın del Rey, and G. Rodr ́ıguez S ́anchez. Modeling epidemics using cellular automata. Applied Mathematics and Computation, 186(1):193 – 202, 2007.
