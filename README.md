### VPNL/invPRF
Code and data accompanying **Poltoratski, Kay, Finzi, and Grill-Spector (2020). Holistic face recognition is an emergent phenomenon of spatial processing in face-selective regions.**

**Abstract:** Spatial processing by receptive fields is a core property of the visual system. However, it is unknown how spatial processing in high-level regions contributes to recognition behavior. As face inversion is thought to disrupt typical ‘holistic’ processing of information in faces, we mapped population receptive fields (pRFs) with upright and inverted faces in the human visual system. In face-selective regions, but not primary visual cortex, pRFs and overall visual field coverage were smaller and shifted downward in response to face inversion. From these measurements, we successfully predicted the relative behavioral detriment of face inversion at different positions in the visual field. This correspondence between neural measurements and behavior demonstrates how spatial processing in face-selective regions may enable holistic perception. These results not only show that spatial processing in high-level visual regions is dynamically used towards recognition, but also suggest a powerful approach for bridging neural computations by receptive fields to behavior. 

_________________________

Contact: Lead author Sonia Poltoratski (sonia09@stanford.edu) and senior author Kalanit Grill-Spector (kalanit@stanford.edu) will handle additional info requests.
_________________________

**Dependencies:** https://github.com/soniapolt/utils + https://github.com/kendrickkay/knkutils.   
Scripts were developed and validated using MATLAB 2019b (Mac OS/Linux). Some require Psychotoolbox 3.0.14 (http://psychtoolbox.org/)   

**Organization:** 

- `/scan-experiment`: stimulus code for the upright/inverted pRF mapping experiment + RSVP letter task
   
   - `/run-output`: .mat files generated during scan runs

- `/scan-analysis`: analysis and plotting code for scan experiment  

    - `/cssFit`: a sample local (non-cluster-compute) version of CSS fitting  

    - `/simulations`: code implementing additional simulations described in figures S4 (noiseSim_), S12 (imOffset_), and S13 (prfRec_)

    - `/coverage`: pre-computed bootstrapped coverage metrics used in coverage density plotting  

    - `/prfSets`: pre-computed and trimmed pRF model fits used in core analyses, as well as alternative stimulus coding fits (figure S12) 
    
    - `/voxPlots`: sample voxelwise results plots, following conventions of figure 1 

- `/scan-eyetracking`: data and analysis code for eye tracking data collected during pRF mapping scans 

- `/behav-experiment`: stimulus and development code for the behavioral FIE experiment, including initial optimal-position simulation (figure 6)     

    - `/data`:	pre-processed eyetracking data from behavioral experiment (converted from eyelink .asc, filtered for fixation breakage)
    
    - `/initialSimulation`:	simulation code and results deriving "optimal" position from initial N=6 data
    
    - `/stims`: stimuli for practice and main experiment	

- `/behav-analysis`: analysis and plotting code for behavioral experiment  

Code underlying each figure is marked via suffix **(e.g. figS1_*)**
Additional code not directly linked to the manuscript: https://github.com/soniapolt/facePRF
