% practicing and running prfRec pilot experiment
% SP 7/2019

s = 'SP'; % set here at the beginning of the session

%%% 1: CD to experiment directory and load images
cd('~/Experiments/prfRec');

if ~exist('alphabet.m','file') addpath(genpath('~/Experiments/Sonia/matlab/utils')); end
load('stims/faceViews.mat'); % now we can do this in the background before starting the experiment

%%% 2: Provide subject with consent form

%%% 3: Show keynote explaining 2-back task, eyetracker calibration.

%%% 4: Run practice trials & eyetracking demo
prfRec2_practice(s,face);

%%% 5: Run 2-back task experiment                              
  % now, you can check params
    [params,expt,block,perf,eyeInit] = prfRec2(s,face);
                                                                                                         
%%% 6: Provide subject with payment form & payment 
%%% ($5 for 30 mins, $10 for 60 mins)
%%% 6b) Log payment in XLS sheet

%%% 7: Optional: Push data to server.
% 7a) Turn on Wi-fi
% 7b) In finder, press COMMAND + K and Connect to: smb://oak-smb-fs-kalanit.stanford.edu/group/biac2/kgs/projects
% 7c) Run pushData commands
pushData('~/Experiments/prfRec/data','/Volumes/projects/behavFIE/prfRec/data','.mat');
pushData('~/Experiments/prfRec/edfs','/Volumes/projects/behavFIE/prfRec/edfs','.edf');

%%% 8: Collect subject demographics
logName = '/Volumes/projects/behavFIE/prfRec/prfRec2_log.mat';
addToLog(logName);

%%% 9: Optional: Make plots
% figure;prfRec_plot([perf.percCorrect],'Percent Correct',params.locNames,[4 5 6 1 2 3]);
% figure;prfRec_plot([perf.dprime],'d Prime',params.locNames,[4 5 6 1 2 3]);
