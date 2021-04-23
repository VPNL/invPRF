% loads & plots distributions of XY changes/anything else per each subject

clear all; %close all;

subjs = prfSubjs;
expt = 'fixPRF';

convertDVA = 1;
flipX = 1;

minR2 = 'r2-20';          % cutoff for vox selection
ROIs= ['V1' standardROIs('face')];

whichStim = 'outline';
whichModel = 'cssExp1';%'kayCSS'; % kayCSS is main text, cssExp1 is non-CSS pRF model
whichM = 'mean';

plotPars = {'X' 'Y' 'gain' 'size' 'r2'};
parTitles = {'X Estim.' 'Y Estim.' 'Gain Estim' 'Size [Sigma/sqrt(N)] (dva)' 'r2'};

hems = {'lh' 'rh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
nBins = 20; % histogram bins

% now we load in the data from both hemispheres, and threshold across
prfSet = ['prfSets/fixPRF_' whichModel '_outline_' hemText(hems) '_r2-20.mat']; % pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems);
load(prfSet);

niceFig([.1 .4 .2*length(plotPars) .6]);

for p = 1:length(plotPars)
    plotPar = plotPars{p}; parTitle = parTitles{p};
    subplot(1,length(plotPars),p);
    
    allPars = nan(2*length(ROIs),length(subjs)); n=1;
    colors = [];
    for r = 1:length(ROIs)
        
        ROInum = cellNum(ROIs{r},info.ROIs);
        
        if strcmp(info.expt,'fixPRF') cOrder = [2 1]; else cOrder = 1:length(subj(1).roi(1).fits); end
        for c = cOrder
            for s = 1:length(subjs)
                fits = subj(s).roi(ROInum).fits(c);
                try
                    eval(['allPars(n,s) = nan' whichM '(getPar(plotPar,fits,1,flipX));']);
                catch allPars(n,s) = NaN; end % for missing values
            end
            allFactors{n} = [ROIs{r} '-' fits.cond(1:3)];
            n=n+1;
            if c == 1 mult = .5; else mult = 1; end
            colors = [roiColors(ROIs{r})*mult; colors];
        end
    end
    
    if containsTxt(plotPar,'gain') cutY = [0 6];
    elseif containsTxt(plotPar,'r2') cutY=[0 100];
    else cutY = []; end
    set(gca,'TickDir','out'); 
    niceBoxPlusGrouped(allPars',ROIs,{roi(1).fits(cOrder).cond},flipud(colors),cutY,1,0);
    title({plotPar;strTogether({roi(1).fits(cOrder).cond})}); titleText = [expt ' (' hemText(hems) '), voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
    axis square;
    superTitle(titleText,titleSize,.025);
end

if onLaptop playSound; end