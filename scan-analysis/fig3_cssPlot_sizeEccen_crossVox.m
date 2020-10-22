% loads & plots model properties (like size x eccen) along the visual heirarchy for cssFit properties
% 7/27/18: change from lsline to fitline2derror
% 3/10/20: change from fitline2derror/scatterline to bootstrapping over
% fitl1line

clear all; close all;


expt = 'fixPRF';

minR2 = 50;          % cutoff for vox selection
ROIs= ['V1' standardROIs('face')];
scatterVox = 300;
doAlpha = 0;
bootMethod = '1l';%'scatterline';%'1l'

whichStim = 'outline';
whichModel = 'kayCSS';

hems = {'rh' 'lh'};
yl = 6;

boot.numIter = 100;
boot.sampleVox = .8;
boot.binThresh = 10; % need n+1 voxels to bootstrap errorbar
boot.fitRange = [.25:.25:6];%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;

load(['prfSets/fixPRF_kayCSS_outline_' hemText(hems) '_r2-50.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f(1) = niceFig([.1 .1 .8 .4],fontSize,1);
%f(2) = niceFig([.1 .1 .7 .4],fontSize,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: vox scatter + fit lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = 1:length(ROIs)
    
    ROInum = cellNum(ROIs{r},info.ROIs);
    fits = roi(ROInum).fits;
    
    if scatterVox < 1 sv = round(length(fits(1).vox)*scatterVox);
    elseif scatterVox > length(fits(1).vox) sv = length(fits(1).vox);
    else sv = scatterVox; end
    sv = datasample([1:length(fits(1).vox)],sv,'Replace',false);
    
    for c = 1:length(fits)
        figure(f(1)); subplot(1,length(ROIs),r);
        if c == 1 mult = .15; else mult = .5; end
        if doAlpha
            alpha = mat2gray([fits(c).vox(sv).r2]);
        else alpha = ones(1,length(sv));end
        hold on;
        s{c} = scatterAlpha([fits(c).vox(sv).eccen],[fits(c).vox(sv).size],alpha,white*mult,2); hold on;
        
        
        %%% bootstrapping stage
        switch bootMethod
            case '1l'
                
                h = NaN(boot.numIter,2); R2 = NaN(boot.numIter,1);
                
                % points at which we'll compute the bootstrapped std
                [N,edges] = histcounts([fits(c).vox.eccen],boot.fitRange);
                fitRange = [boot.fitRange(N>boot.binThresh)' ...
                    ones(length(boot.fitRange(N>boot.binThresh)),1)]; fitOut = NaN(boot.numIter,length(fitRange));
                
                %tic
                parfor b = 1:boot.numIter
                    v = datasample([1:length(fits(1).vox)],round(length(fits(1).vox)*boot.sampleVox),'Replace',true);
                    x = [fits(c).vox(v).eccen]';
                    y = [fits(c).vox(v).size]';
                    X = [x ones(length(v),1)];
                    [h(b,:),R2(b)] = fitl1line(X,y);
                    
                    fitOut(b,:) = fitRange*h(b,:)';
                end
                %toc
                
                % full line fit
                x = [fits(c).vox.eccen]';
                y = [fits(c).vox.size]';
                X = [x ones(length(fits(c).vox),1)];
                
                [hFull,R2Full] = fitl1line(X,y);
                h1 = plot(boot.fitRange,[boot.fitRange' ones(length(boot.fitRange),1)]*hFull','Color',roiColors(ROIs{r})*c*.5); hold on;
                [ci,med] = CI(fitOut);
                
                % construct errorbar
                plotErr = NaN(2,length(boot.fitRange));
                plotErr(:,find(N>boot.binThresh)) = ci;
                
                hold on; e = errorbar3(boot.fitRange,[[boot.fitRange' ones(length(boot.fitRange),1)]*hFull']',plotErr,'v',roiColors(ROIs{r})*c*.5); set(e,'FaceAlpha',1);
                
            case 'scatterline'
                x = [fits(c).vox.eccen]';
                y = [fits(c).vox.size]';
                %X = [x ones(length(fits(c).vox),1)];
                
                
                [N,edges] = histcounts(y,boot.fitRange);
                [errorObj,lineObj,mn,se] = scatterline(x,y,boot.fitRange(N>boot.binThresh),NaN,1000,roiColors(ROIs{r})*mult,2,1); hold on;
                extX = [.25:.1:6];
                extLine = plot(extX, extendLine(lineObj,extX),'Color',roiColors(ROIs{r})*mult,'LineWidth',.5,'LineStyle','--');
                set(lineObj,'LineWidth',.5); set(extLine,'LineStyle',':');
                set(errorObj,'FaceAlpha',.2);
                
        end
        
        
        if c == 2
            xl = xlim; xlim([0.25 6]); y = ylim; if containsTxt(ROIs{r},'faces') ylim([.25 yl]); else ylim([.25 6]);end
            set(gca,'TickDir','out');
            xlabel('Eccen (dva)','FontSize',fontSize); ylabel(['Size (2*SD/sqrt(N)) (dva)'],'FontSize',fontSize);
            axis square;
            title(ROIs{r});
        end
        
        %superTitle(titleText,titleSize);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % figure 2: fit lines for all ROIs (separated by face/evc)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
%         figure(f(2));
%         %if containsTxt(ROIs{r},'faces') subplot(1,2,2); else subplot(1,2,1); end % summary figure across ROIs
%         subplot(1,2,c);
%         switch bootMethod
%             case '1l'
%                 hold on; h1 = plot(boot.fitRange,[boot.fitRange' ones(length(boot.fitRange),1)]*hFull','Color',roiColors(ROIs{r})*c*.5);
%                 hold on; lo{r} = errorbar3(boot.fitRange,[[boot.fitRange' ones(length(boot.fitRange),1)]*hFull']',plotErr,'v',roiColors(ROIs{r})*c*.5); set(lo{r},'FaceAlpha',1);
%             case 'scatterline'
%                 hold on; [eo{r},lo{r},mn,se] = scatterline([fits(c).vox.eccen]',[fits(c).vox.eccen]',boot.fitRange(N>boot.binThresh),NaN,1000,roiColors(ROIs{r})*mult,2,1); hold on;
%                 extLine = plot(extX, extendLine(lo{r},extX),'Color',roiColors(ROIs{r})*mult,'LineWidth',.5);
%                 set(lo{r},'LineWidth',.5); set(extLine,'LineStyle',':');
%                 set(eo{r},'FaceAlpha',.2);
%         end
    end
end


playSound;