% feb 2020 version that extends noiseSim_iter2.m to manuscript-version
% simulation
% takes X random voxels from a specified ROI across all subjects, scales
% betas and iteratively adds noise to match upright R2 to inverted R2,
% compares resulting estimates to 

clear all; close all;

rois = {'mFus_faces'};
numDraws = 50;

for rr = 1:length(rois)
 
% experiment and session
sim.expt = 'fixPRF';
sim.numVox = .8;
sim.whichModel = 'kayCSS'; sim.whichStim = 'outline'; sim.minR2 = 20;
sim.simSuffix = 'scaleIterNoise';
sim.iterStep = .9; % starting value for noise iteration
sim.r2thresh = 1;   % in r2 units, difference between conditions that we'll allow

sim.ROI = rois{rr};%'pFus_faces';
sim.hem = 'rh';
tic

% retrieve session
sim.baseCond = 2; sim.compCond = 1;

% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),sim.expt,sim.minR2,sim.whichStim,sim.whichModel,{sim.hem}));
ROInum = cellNum(sim.ROI,info.ROIs);

fits = roi(ROInum).fits; % since we're looking at one ROI at a time here, simplify
if sim.numVox < 1 sim.numVox = ceil(sim.numVox*length(fits(1).vox)); end % if we've given a proportion of voxels

 for d =  2:numDraws
sim.voxInd = datasample(1:length(fits(1).vox),sim.numVox,'Replace',false); % indices of voxels from initial pRFset
sim.base = fits(sim.baseCond).vox(sim.voxInd);
sim.comp = fits(sim.compCond).vox(sim.voxInd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sample vox, add noise, refit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%sim.vox = struct;

for v = 1:sim.numVox
    % which session does this voxel come from?
    [~,stimN]=fileparts(sim.base(v).stim);
    sim.vox(v).session = stimN(1:8);
    
    % load data & stim - upright
    load([dirOf(pwd) 'cssFit/' sim.base(v).stim]); % downsampled and reshaped
    stimulus = reshape(condAvg(:,:,fits(sim.baseCond).condNums),im.size*im.size,length(fits(sim.baseCond).condNums))';
    
    % the parameters of the CSS model are [R C S G N] where
    %   R is the row index of the center of the 2D Gaussian
    %   C is the column index of the center of the 2D Gaussian
    %   S is the standard deviation of the 2D Gaussian
    %   G is a gain parameter
    %   N is the exponent of the power-law nonlinearity
    %   B is the baseline shift
    
    switch sim.whichModel
        case 'cssExpN'
            sim.exptN = .2;
            sim.modelName = ['singleGauss CSS + fixed exponent: ' num2str(sim.exptN)];
            [fits.parNames] = deal({'Y','X','sd','gain'});
            [fits.sim.exptN] = deal(sim.exptN);
            [modelfun, model, metric, resampling] = init_cssexptN(sim.exptN,im.size,sim.hem,im.ppd);
 
        case 'kayCSS'
            sim.modelName = 'kay CSS, no baseline shift';
            [fits.parNames] = deal({'Y','X','sd','gain','exp'});
            [modelfun, model, metric, resampling] = init_kayCSS(im.size,sim.hem,im.ppd);
    end
    
    if containsTxt(sim.simSuffix,'scale') % are we rescaling the beta amps?
        sim.vox(v).scale = range(sim.base(v).betas)/range(sim.comp(v).betas);
    else % oh no?
        sim.vox(v).scale =1;
    end
    
    sim.vox(v).scaleData = sim.base(v).betas./sim.vox(v).scale;
    
    sim.vox(v).noiseLevel = abs(sim.base(v).r2-sim.comp(v).r2)/100;
    
    while 1 % add noise iteratively until r2 is ~ matched
        i = 1;
        fprintf('Now simulating with noise level = %.4f...\n',sim.vox(v).noiseLevel);
        sim.vox(v).addNoise = sim.vox(v).noiseLevel*randn(length(sim.base(v).betas),1);
        sim.vox(v).simData = sim.vox(v).scaleData+sim.vox(v).addNoise;
        
        opt = struct( ...
            'stimulus',    stimulus, ...
            'data',        sim.vox(v).simData, ...
            'model',       {model}, ...
            'resampling',  resampling, ...
            'metric',      metric, ...
            'dosave',      'modelfit',...
            'optimoptions', {{'Display' 'off' 'UseParallel' 0}}...
            );
        
        %%% fit the model
        results = fitnonlinearmodel(opt);
        sim.vox(v).results = results;
        sim.vox(v).params = results.params;
        sim.vox(v).r2 = results.trainperformance;
        sim.vox(v).modelfit = results.modelfit;
        sim.vox(v).iter = 1; % success
        i=i+1;
    if abs(sim.vox(v).r2-sim.comp(v).r2)<=sim.r2thresh
        fprintf('Voxel %d finished fitting in %d iterations...\n',v,i);
        break
    elseif sim.vox(v).noiseLevel < .001 % something weird happened (likely comp cond was less noisy to start with)
        sim.vox(v).iter = 0; % not success
        fprintf('Voxel %d could not fit.\n',v);
        break
    elseif sim.vox(v).r2 < sim.comp(v).r2 % less noise
            sim.vox(v).noiseLevel =   sim.vox(v).noiseLevel * sim.iterStep;
        else sim.vox(v).noiseLevel =   sim.vox(v).noiseLevel/sim.iterStep; % more noise
        end
    end
   
end

save(['iterNoiseSim/' sim.simSuffix '_' sim.expt '_' sim.whichModel '_' sim.hem sim.ROI '_' num2str(d) '.mat'],'sim');
timeFit = toc;
fprintf('%s modelFit finished %s%s at %s, took %s mins\n',sim.whichModel,sim.hem,sim.ROI,datestr(now),num2str(timeFit/60));

end
end