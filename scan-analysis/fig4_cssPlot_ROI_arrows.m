% plots XY shift for ROIs, both bilaterally and for each individual hem

clear all; close all;

subjs = prfSubjs;
expt = 'fixPRF';

minR2 = 'r2-50';        % cutoff for vox selection
ROIs= {'mFus_faces'};% standardROIs('face')];%standardROIs;%['hV4' standardROIs('face')]; %{'ML' 'PL'};%{};%('face')
% manual set of baseCond + compConds
% [baseCond, compCond], more flexibly defined
base = 2; comp = 1;

saveFig = 0;

sampleVox = 300; % number of voxels to grab for plotting

whichStim = 'outline';%'edge';%'binary';%'internal';%
whichModel = 'kayCSS';%'inflipCSSn';%'kayCSS';%'tempCSSn';%'%cssShift';%
fitSuffix = '';

hems = {'lh' 'rh'};%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fontSize = 12; titleSize = 14;
niceFig([.1 .1 .8 .8],12,1);numPlots = [3,length(ROIs)*2];
for r = 1:length(ROIs)
    for h = 1:2
        %load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems(h),fitSuffix));
        load(['prfSets/fixPRF_kayCSS_outline_' hems{h} '_r2-50.mat']);
        
        ROInum = cellNum(ROIs{r},info.ROIs);
        fits = roi(ROInum).fits;
        
        if sampleVox < 1 sv = round(length(fits(1).vox)*sampleVox);
        elseif sampleVox > length(fits(1).vox) sv = length(fits(1).vox);
        else sv = sampleVox; end
            randVox = datasample([1:length(fits(1).vox)],sv,'Replace',false);
            for c = 1:length(fits) fits(c).vox = fits(c).vox(randVox); end % we do this so we're grabbing the same randomized voxels in hem and bilat plots
            if h == 1 bFits = fits; else for c = 1:length(fits) bFits(c).vox = [bFits(c).vox fits(c).vox]; end % assemble bilateral fits
            end
            
            subplot(numPlots(1),numPlots(2),2*(r-1)+h);
            
            plotXYshift(fits(base).vox,fits(comp).vox,fits(1).ppd,fits(1).res,3);
            drawnow;
            
            ylabel([fits(base).cond ' to ' fits(comp).cond],'FontSize',12,'FontWeight','bold');
            set(get(gca,'YLabel'),'Visible','on');
            title([hems{h} ROIs{r}]);
            
            if h == 2
                subplot(numPlots(1),numPlots(2),[2*(r-1)+h-1+numPlots(2) 2*(r-1)+h+numPlots(2) ...
                    2*(r-1)+h-1+2*numPlots(2) 2*(r-1)+h+2*numPlots(2)]);
                
                plotXYshift(bFits(base).vox,bFits(comp).vox,fits(1).ppd,fits(1).res);
                drawnow;
                
                
                % [left bottom width height]
                pos = get(gca,'Position');
                legSz = .08;
                axes('Position',[pos(1) pos(2)+pos(4)-2*legSz legSz legSz]);
                drawcolorbarcircular(cmapang,1);
                title([ROIs{r} ' (' num2str(2*length(randVox)) 'total voxels R^2 > ' num2str(minR2) ')']);
                
            end
        end
    end
    
%     if saveFig
%         txt = [whichModel '_' whichStim '_' minR2 '_voxSample' num2str(sampleVox) '_XYshift'];
%         niceSave([dirOf(pwd) 'figures/' expt '/arrows/' ],txt,[],[],{'png' 'svg'}); % just save pngs, since these can be generated pretty quickly
%     end
    
    superTitle([whichModel '-' whichStim ' , Subjs: ' strTogether(info.subjs) ...
        ' (' num2str(sampleVox) 'rand-sampled voxels R^2 > ' num2str(minR2) ')'],12,.97);
    
    if onLaptop playSound; end