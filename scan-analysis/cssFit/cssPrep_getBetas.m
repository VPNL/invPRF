function [outputFile] = cssPrep_getBetas(expt,thisROI,session,numScans,dt,recomp)
%clear all;
%thisROI = vpnlROI('lh_V1'); expt = 'fixPRF'; session = 'AS180807'; numScans = 9;

if ~exist('dt','var')|| isempty(dt);
    dt = 4; % default - slice time corrected, motion corrected, single-session.
end

initDir = pwd;
mainDir = fullfile('/share','kalanit','biac2','kgs', 'projects', 'invPRF'); % path to raw data
whichBetas = 'GLM'; % 'GLM' or 'deconv'

scan1 = 1; % first scan in this group

outputDir = fullfile(mainDir,'cssFit','data',expt,session);
checkDir(outputDir);

output = struct('ROIname',thisROI,'betas',[],'sems',[]);

cd( fullfile(mainDir,expt, session) );
fprintf('*** Session %s *** \n', session)

%enforce consistent preprocessing / event-related parameters
if ~exist([mainDir expt '/' session '/Inplane/ROIs/' thisROI '.mat'])
    hI = initHiddenInplane(dt, scan1);
    hI = roiLoadVol2Inplane(hI,['3DAnatomy/ROIs/' thisROI],1);
end

hI = initHiddenInplane(dt, scan1,thisROI);

% groups scans & assigns parfiles, if we haven't already
glm_groupScans(dt,session,expt);

params = er_getParams(hI);

if strcmp(whichBetas,'GLM') ==1
    params.glmHRF = 3; % flag for spm HIRF
    params.ampType = 'betas';
else params.ampType = 'deconvolved'; end

params.detrend = 1; % high-pass filter (low-pass is blurring, below)
params.inhomoCorrect= 1; % 1 = divide by mean, 2 = subtract null
params.temporalNormalization = 0; % no matching 1st temporal frame
params.eventsPerBlock = 2; % TRs per stim block
params.framePeriod = 2; % fixes a bug introduced in the pre-processing code for fixPRF subjs (EM, AS, TH, MN)
er_setParams(hI, params);

if ~exist('recomp','var')
    setpref('VISTA', 'recomputeVoxData',0);
else setpref('VISTA', 'recomputeVoxData',1); end

% this grabs the data for each voxel in a 'multivoxel' structure
mv = mv_init(hI);

% run glm at the level of each voxel
mv=mv_applyGlm(mv); %mv.glm.betas is what you are going to want to use

% remove NaN voxels
remove = [];
for n= 1:length(mv.glm.varianceExplained)
    if isnan(mv.glm.varianceExplained(n))
        fprintf('Removing voxel #%d/%d for NaN values...\n',n,length(mv.glm.varianceExplained))
        remove(end+1) = n;
    end
end

% organize GLM output, finish removing
output.betas = squeeze(mv.glm.betas(:,1:end-numScans,:));
output.sems = squeeze(mv.glm.sems(:,1:end-numScans,:));
if length(remove)>0
    output.betas(:,remove) = [];
    output.sems(:,remove) = [];
end

outputFile = [outputDir '/' session '_' thisROI '.mat'];
save(outputFile,'mv','output');
setpref('VISTA', 'recomputeVoxData',0);
cd(initDir);
%end
