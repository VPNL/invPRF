% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs = prfSubjs;%{'MG' 'JG' 'TH' 'EM' 'DF' 'SP'};%
expt = 'fixPRF';

saveFig =1;
convertDVA = 1;

minR2 = 'r2-20';%['perc-50'];          % cutoff for vox selection
ROIs= ['V1' standardROIs('face')];%{'mFus_faces'};%standardROIs;%('face+')

whichStim = 'outline';%'photo';%'eyes';%'internal';%
whichModel = 'kayCSS';%'flipCSSn';%'cssExpN';%'cssShift';%
whichM = 3; % 1 = mean, 2 = mode/peak, 3 = median
plotPars = {'r2'};%{'Y' 'size'};%{'gain' 'r2' 'Y' 'X' 'size' };%{'Y'};%
parTitles = {'Estimated R^{2}'};%{'Y' 'size [sigma/sqrt(N)] (dva)'}%{ 'Gain Estim' 'Estimated R^{2}' 'Y Estim.' 'X Estim.' 'Size [2xSD/sqrt(N)] (dva)'};%{'Y estim'};%
plotType = {'scatter' 'scatter'};%{'box' 'scatter' 'distr' 'distr' 'scatter' };%{'distr'};% % 1 = boxplot, 2 = distr, 3 = scatter

hems = {'rh' 'lh'};
fitSuffix = '';%'_orig';%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
nBins = 20; % histogram bins

% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

if length(subjNum) == 1 roi = subj(subjNum).roi; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: distribution of parameters for this ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:length(plotPars)
    
titleText = [whichModel ' ' parTitles{p} ', Subj: '];
titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if saveFig && onLaptop figSize = [0.1 0.1 .6 .4]; else figSize = [.2 .1 8 .8]; end
    niceFig(figSize,fontSize,1);
    numPlots = [1 ceil(length(ROIs)/1)];pl = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for r = 1:length(ROIs)
    fits = roi(ROInum(r)).fits;
    subplot(numPlots(1),numPlots(2),pl)
    
    for c = 1:length(roi(1).fits)
        
        cPars{c} = getPar(plotPars{p},fits(c),1);
        if containsTxt(plotPars{p},'size')
                cPars{c}= cPars{c}./2;    
                end
%         parNum = cellNum(plotPars{p},fits(1).parNames);
%             if ~isempty(parNum)
%                 pars = vertcat(fits(c).vox.params);
%                 cPars{c} = pars(:,parNum)';
%             else
%                 eval(['cPars{c} = [fits(c).vox.' plotPars{p} '];']);  end
%         
%          
%         % rescale some parameters so that they are in DVA units and
%         % centered around zero (center of screen)
%         if convertDVA && containsTxt(plotPars{p},'Y') || containsTxt(plotPars{p},'X') || containsTxt(plotPars{p},'sd')
%         if ~containsTxt(plotPars{p},'sd') % don't re-center the SD
%             cPars{c} = fits(1).res-cPars{c}-roi(1).fits(1).res/2;
%         end
%         cPars{c} = cPars{c}./roi(1).fits(1).ppd;
%         end
    end 
    
    switch plotType{p}
            case 'box'
                if containsTxt(plotPars{p},'gain') 
                    niceBoxplot([cPars{1};cPars{2}]',{fits(1).cond fits(2).cond},1,[condColors(4);condColors(2)],[0 5]);
                   % ylim([0 5]);
                elseif containsTxt(plotPars{p},'r2')                     
                    niceBoxplot([cPars{1};cPars{2}]',{fits(1).cond fits(2).cond},1,[condColors(4);condColors(2)],[minR2 100]);
                else
                niceBoxplot([cPars{1};cPars{2}]',{fits(1).cond fits(2).cond},1,[condColors(4);condColors(2)]);
                end
                if containsTxt(plotPars{p},'size') 
                   ylim([0 5]); 
                end
            case 'distr'
                nBins = 20;
                plotDistr(cPars,1,{fits(1).cond fits(2).cond},nBins,3,1);
                if containsTxt(plotPars{p},'size') 
                   xlim([0 10]); 
                else xlim([-5 5]); end
                
                [pv, od, effectsize] = permutationTest(cPars{1}, cPars{2}, 1000, ...
                     'plotresult', 0, 'showprogress', 0);
                 fprintf('%s Permutation Test, %s: p = %.6f, observed difference = %.2f\n ',...
                     ROIs{r}, plotPars{p},pv,od);
            case 'scatter'
                hold on;
                colors = {roiColors(ROIs{r}); roiColors(ROIs{r})*.5};
                scatterCent2(cPars{1},cPars{2},colors,...
                    fits(1).cond,fits(2).cond,[],fontSize,0,1);
                %hold on; l =lsline; l=fliplr(l);
                if strcmp(plotPars{p},'r2') xlim([20 100]); ylim([20 100]); end
                hold on; xl = get(gca,'xlim');
                hold on; plot(xl,xl,'k:');
        end
        %title([whichModels{t}  '-' allStims{t} ' stim' ' ' plotPars{p}],'fontSize',titleSize,'interpreter','none','FontWeight','bold'); 
        pl = pl+1;
        axis square;
    
    
    title({ROIs{r};parTitles{p}},'fontSize',titleSize,'interpreter','none','FontWeight','bold'); 
    %xlabel(roi(1).fits(1).parNames{p},'fontSize',titleSize);
    end
    superTitle(titleText,titleSize,.025);

if saveFig == 1
    if length(subjs) == 1 txt = ['subj' subjs{1}]; else txt = ['groupN' num2str(length(subjs))]; end
    txt = [plotPars{p} '_' txt ];
    
    if length(hems) == 1
        txt = [txt '_' hems{1}]; end
        txt = [plotType{p} '_' whichModel '_' whichStim '_' txt];
    niceSave([dirOf(pwd) 'figures/' expt '/crossVox/'],txt,[],[],{'png' 'svg'}); % just save pngs, since these can be generated pretty quickly
end
end
if onLaptop playSound; end