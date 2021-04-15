% loads & plots distributions of XY changes/anything else per each subject
%%%% this is the new supplement fig for the revision

clear all; close all;

subjs = prfSubjs;
expt = 'fixPRF';

saveFig =1;
convertDVA = 1;

minR2 = 'r2-20';          % cutoff for vox selection
ROIs= ['V1' standardROIs('face')]

whichStim = 'outline';%
whichModel = 'kayCSS';%
whichM = 3; % 1 = mean, 2 = mode/peak, 3 = median
mS = {'mean' 'mode' 'median'};
plotPars = {'gain' 'r2' 'Y' 'X' 'size'};%
parTitles = {'Gain Estim' 'Estimated R^{2}' 'Y Estim.' 'X Estim.' 'Size [2xSD/sqrt(N)] (dva)'};

hems = {'rh' 'lh'};
fitSuffix = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
nBins = 20; % histogram bins

% now we load in the data from both hemispheres, and threshold across
load(['prfSets/fixPRF_kayCSS_outline_' hemText(hems) '_r2-20.mat']);
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

if length(subjNum) == 1 roi = subj(subjNum).roi; else subj = subj(subjNum); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: distribution of parameters for this ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:length(plotPars)
    
    titleText = [whichModel ' ' parTitles{p} ', Subj: '];
    titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if saveFig && onLaptop figSize = [.1 .1 .8 .4]; else figSize = [.2 .1 8 .8]; end
    niceFig(figSize,fontSize);
    numPlots = [1 ceil(length(ROIs))];pl = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for r = 1:length(ROIs)
        
        subplot(numPlots(1),numPlots(2),pl)
        
        sPars = []; mPars = [];
        for c = 1:length(roi(1).fits)
            for s = 1:length(subjNum)
                fits = subj(subjNum(s)).roi(ROInum(r)).fits(c);
                
                try
                cPars{c} = getPar(plotPars{p},fits,1);
                catch cPars{c} = NaN; end % missing ROIs
                
                % aggregate across subjects
                sPars{s,c} = cPars{c}; % full distribution for this subject
                eval(['mPars(s,c) = nan' mS{whichM} '(cPars{c});']); % mean value for this subject
            end
        end
        hues = [1 .5];
        if containsTxt(plotPars{p},'size')
            xl = [0 5] ;
        else xl = [-2 2]; end
        for n = 1:length(sPars)
            dPars{n} = sPars{n,1}-sPars{n,2};
        end
        [h,meds, normcounts] = plotMeanDistr(dPars,nBins,roiColors(ROIs(r)),1); %xlim(xl);
        set(gca,'TickDir','out');
        pl = pl+1;
        axis square;
        
        [h,pv,ci] = ttest(mPars(:,2)-mPars(:,1));
        title({ROIs{r};parTitles{p};['ttest on means: p=' num2str(pv)]},'fontSize',titleSize,'interpreter','none','FontWeight','bold');
    end
    superTitle(titleText,titleSize,.025);
    
%     if saveFig == 1
%         if length(subjs) == 1 txt = ['subj' subjs{1}]; else txt = ['groupN' num2str(length(subjs))]; end
%         txt = [plotPars{p} '_' hemText(hems) '_' txt];
%         niceSave([dirOf(pwd) 'figures/' expt '/deltaDist/'],txt,[],[],{'png' 'svg'}); % just save pngs, since these can be generated pretty quickly
%     end
end
if onLaptop playSound; end