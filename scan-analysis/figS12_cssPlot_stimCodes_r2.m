% (fairly hack-y) viz of the R2 for all ROIs across stim codings

clear all; close all;

subjs = prfSubjs;
expt = 'fixPRF';

saveFig = 0;
convertDVA = 1;

minR2 = 'r2-0';
ROIs= {'V1' 'mFus_faces'};%standardROIs('face')];

modelStims = {'disk','outline','contrast', 'internal','eyes'};
whichM = 3; % 1 = mean, 2 = mode/peak, 3 = median
plotPars = {'r2'};
plotType = {'box'};

hems = {'rh' 'lh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 8; titleSize = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: distribution of parameters for this ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = 1:length(ROIs)
    niceFig('wide',fontSize,1);
    numPlots = [1 length(modelStims)];pl = 1;
    
    p=1; % plot just one par
    for t = 1:length(modelStims)
        whichStim = modelStims{t};
        load(['prfSets/fixPRF_kayCSS_' whichStim '_bilat_r2-50.mat']);
        ROInum = cellNum(ROIs,info.ROIs);
        subjNum = cellNum(subjs,info.subjs);
        
        if length(subjNum) == 1 roi = subj(subjNum).roi; end
        
        subplot(numPlots(1),numPlots(2),t)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fits = roi(ROInum(r)).fits;
        %subplot(numPlots(1),numPlots(2),pl)
        
        for c = 1:length(roi(1).fits)
            
            cPars{c} = getPar(plotPars{p},fits(c),1);
        end
        niceBoxplot([cPars{1};cPars{2}]',{fits(1).cond fits(2).cond},1,[condColors(4);condColors(2)]);
        ylim([0 100]);
        
        pl = pl+1;
        axis square;
        
        title(whichStim,'fontSize',titleSize,'interpreter','none','FontWeight','bold');
    end
    superTitle(ROIs{r},10,.05)
    if saveFig == 1
        
        txt = ['r2_modelComp_' ROIs{r}];
        niceSave([dirOf(pwd) 'figures/' expt '/modelComp/'],txt,[],[],{'png' 'svg'}); % just save pngs, since these can be generated pretty quickly
    end
end
if onLaptop playSound; end