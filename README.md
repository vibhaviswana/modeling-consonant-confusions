# Scene-analysis modeling of consonant confusions

This repository contains the code used in 
[Viswanathan et al. (2021)](https://www.biorxiv.org/content/10.1101/2021.09.06.459159v1).

Viswanathan, V., Shinn-Cunningham, B. G., & Heinz, M. G. (2021). Speech
categorization reveals the role of early-stage temporal-coherence processing
in auditory scene analysis. bioRxiv. DOI: https://doi.org/10.1101/2021.09.06.459159

## Basic usage

To use the modeling code in this repository, CV syllable recordings in different masking 
conditions will need to be supplied. The code should be run using the following steps: 

- Run ```ANexperiment.m```. This runs the Bruce et al. (2018) auditory-nerve model on the 
different audio stimuli and saves outputs to ```Bruce2018ANoutput_CVstim.mat```.
- To model within-channel modulation masking, run ```withinChannelModel.m```, followed by
```calibration_and_prediction_withinChannelModel.m```.
- To model temporal-coherence-based across-channel modulation interference, run 
```acrossChannelModel.m```, followed by ```calibration_and_prediction_acrossChannelModel.m```.

---
**Note**

The code calls the function ```modFbank_v3.m```, which is available as part of the [Auditory Modeling Toolbox; AMT](https://amtoolbox.org) (see the AMT model ```joergensen2013.m```). The Bruce et al. (2018) model is also available in AMT.

---
