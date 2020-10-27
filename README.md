### VPNL/invPRF
Code and data accompanying **Poltoratski, Kay, Finzi, and Grill-Spector (2020). Holistic face recognition is an emergent phenomenon of spatial integration in face-selective regions.**

**Abstract:** Spatial processing by receptive fields is a core property of the visual system. However, it is unknown how spatial coding in high-level regions contributes to recognition behavior. As face inversion is thought to disrupt typical ‘holistic’ processing of information in faces, we mapped population receptive fields (pRFs) with upright and inverted faces in the human visual system. In face-selective regions, but not primary visual cortex, pRFs and overall visual field coverage were smaller and shifted downward in response to face inversion. From these measurements, we successfully predicted the relative behavioral detriment of face inversion at different positions in the visual field. This correspondence between neural measurements and behavior demonstrates how spatial integration in face-selective regions enables holistic processing. These results not only show that spatial processing in high-level visual regions is dynamically used towards recognition, but also suggest a powerful approach for bridging neural computations by receptive fields to behavior. 

_________________________

Contact: Lead author Sonia Poltoratski (sonia09@stanford.edu) and senior author Kalanit Grill-Spector (kalanit@stanford.edu) will handle additional info requests.
_________________________

**Dependencies:** https://github.com/soniapolt/utils + https://github.com/kendrickkay/knkutils.   
Scripts were developed and validated using MATLAB 2019b (Mac OS/Linux). Some require Psychotoolbox 3.0.14 (http://psychtoolbox.org/)   

**Organization:** 

- `/scan-experiment`: stimulus code for the upright/inverted pRF mapping experiment + RSVP letter task 

- `/scan-analysis`: analysis and plotting code for scan experiment  

    - `/cssFit`: a sample local (non-cluster-compute) version of CSS fitting  

    - `/simulations`: code implementing additional simulations described in figures S2 (noiseSim_), S7 (prfRec_), and S8 (imOffset_)  

    - `/coverage`: pre-computed bootstrapped coverage metrics used in coverage density plotting  

    - `/prfSets`: pre-computed and trimmed pRF model fits used in core analyses  

- `/behav-experiment`: stimulus and development code for the behavioral FIE experiment, including initial optimal-position simulation (figure 6)     

    - `/data`:	pre-processed eyetracking data (converted from eyelink .asc, filtered for fixation breakage)  

- `/behav-analysis`: analysis and plotting code for behavioral experiment  

Code underlying each figure is marked via suffix **(e.g. figS1_*)**
Additional code not directly linked to the manuscript: https://github.com/soniapolt/facePRF
