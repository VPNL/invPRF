% analysis script for initial import and reading of eyelink files
% will produce plots of eye position for each trial, condition-wise
% averages
% updated for prfRec2

%%% ET NAME: always in the format expt_subj_session - will be parsed later


%for dc = [1 0]
    
driftCorr = 1; % binary, whether we apply driftCorrection to the resulting plots

subjs = prfSubjs;

for s = 1:length(subjs)

% set dirs & files
exptDir = [dirOf(pwd) 'prfRec']; 
subj = subjs{s}; sess = prfRecSess(subj);%'190822';%
doParse = 1;

% determine other dirs & files
dataName = sprintf('prfRec2_%s',subj);
etName = sprintf('prfRec2_%s_%s',subj,sess);

ascFile = [exptDir '/ascs/' etName '.asc']; 
dataFile = [exptDir '/data/' dataName  '.mat']; 
matFile = [exptDir '/mat/' etName '.mat'];
figDir =  [dirOf(exptDir) 'analysis/figures/' etName '/'];

% read asc file if necessary - this should be done at the end of the
% eyetracking session
if ~exist(ascFile,'file') edfRead([exptDir '/edfs/' etName '.edf'],[exptDir '/ascs/']); end

% parse ASC file to event & sample info
if ~exist(matFile,'file') || doParse [trial,info] = ascParse(ascFile,dataFile,matFile); else load(matFile); end
info

% plot trials & condition averages
eyetrackPlots(exptDir,etName,dataName,driftCorr,[1 0 1]);

% output some performance - currently taken from orig expt
load(dataFile); 
titleText = sprintf('Session: %s %d, Task: %Stig WM\n\n',subj,sess);
labels = {'LowerLeft' 'Center' 'UpperRight'}; order = [4 5 6 1 2 3];

figure;prfRec_plot([perf.hitRate],'Hit Rate',labels,order); title(titleText); niceSave(figDir,'behav_hitRate');
figure;prfRec_plot([perf.dprime],'d Prime',labels,order);title(titleText);niceSave(figDir,'behav_dPrime');

cd(exptDir);

playSound;


end
%end