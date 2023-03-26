This is a set of computational models associated with the paper "Genetic control of cell layer interactions in plants via tissue mechanics", by Kelly-Bellow et al., to appear in Science, probably in 2023.

There are five models:

GPT\_ArabidopsisTissueCrack\_20230319  GPT\_Utricularia\_20230319  GPT\_UtriculariaSolid\_20230319  GPT\_UtriculariaThirdBladesMissing\_20230319  
80-HypocotylTissueStressSmall

The first four of these model biological tissues as continuous materials subject to elasticity and growth. They can be run using the Matlab package GFtbox, which is available from [GitHub](https://github.com/JIC-Enrico-Coen/GrowthToolbox). GFtbox requires an installation of Matlab, which is commercial software available from [Mathworks](https://Mathworks.com).

The fifth, 80-HypocotylTissueStressSmall, is a simulation of cracking in the hypocotyl based on modelling individual cells. It requires the MorphoDynamX software package, which is available in pre-release form from [morphographx.org](https://morphographx.org/morphodynamx/).