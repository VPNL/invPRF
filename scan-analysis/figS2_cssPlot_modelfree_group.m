% loads & plots some basics re: invPRF fits (css model)
% SP 10/25/17
% _indiv code allows us to plot individual ROIs in multiple plots
% 10/8/18 edit: made more flexible across expriment with different numbers
% of fit conditions

clear all; close all;

% ROIs and subjects
subjs = prfSubjs;
expt = 'fixPRF';
ROIs= ['V1' standardROIs('face')];% %{'PL' 'ML'};%'V1'
sporder = [4 1 3 5 6];

% plotting choices
saveFig = 0;

% fit file
minR2 = 0;          % cutoff for vox selection

whichStim = 'outline';
whichModel = 'kayCSS';
hems = {'rh' 'lh'};

plotOrder = [2 1];
plotX = [1:25; 28:52];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 6; titleSize = 8;

% now we load in the data from both hemispheres, and threshold across
load(['prfSets/fixPRF_kayCSS_outline_' hemText(hems) '_r2-0.mat']);%load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems));
subjs = info.subjs;

niceFig([.1 .1 .8 .8]);

for r = 1:length(ROIs)
    ROI = ROIs{r}
    ROInum = cellNum(ROI,info.ROIs);
    clear betas subjBetas fits condMean
    
    for s = 1:length(subjs)
        
        
        fits = subj(s).roi(ROInum).fits;
        
        
        for c = 1:2
            for v = 1:length(fits(1).vox)
                betas(v,:,:) = reshape(fits(plotOrder(c)).vox(v).betas,sqrt(length(fits(plotOrder(c)).vox(v).betas)),sqrt(length(fits(plotOrder(c)).vox(v).betas)))';
            end
            subjBetas(c,s,:,:) = squeeze(mean(betas));
        end
    end
    
    condMean = squeeze(mean(subjBetas,2));
    
    titleText = ['(N=12) average betas'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % beta value plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%% summary figure
    ax1= subplot(3,2,sporder(r))
    ROI
    cb = [squeeze(condMean(1,:,:))  squeeze(condMean(2,:,:))]
    [min(cb(:)) max(cb(:))]
    plotInSpace([squeeze(condMean(1,:,:))  squeeze(condMean(2,:,:))],ROI,'Upr    -    Inv',1,[min(cb(:)) max(cb(:))]);
    colormap(ax1,parula);
    
    superTitle(titleText,12, .05);
end
if saveFig
    txt = ['allROIs_modelFree-r' num2str(minR2)];
    niceSave([dirOf(pwd) 'figures/' expt '/modelFree/'],txt,[],[],{'png', 'svg'}); % just save pngs, since these can be generated pretty quickly
    
end

if onLaptop playSound; end