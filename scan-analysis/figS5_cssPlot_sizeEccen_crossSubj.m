% loads & plots model properties (like size x eccen) along the visual heirarchy for cssFit properties
% 3/11/20: THIS IS THE FIG 4 CODE

clear all; close all;

expt = 'fixPRF';

minR2 = 20;          % cutoff for vox selection
ROIs= ['V1' standardROIs('face')];
sampleVox = 300;
fitIndivs = 0;
fitRange = [0:.25:8];

whichStim = 'outline';
whichModel = 'kayCSS';%

hems = {'lh' 'rh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
pf = ['prfSets/fixPRF_kayCSS_outline_' hemText(hems) '_r2-20.mat']; % pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems);

load(pf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: vox scatter + fit lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f(2) = niceFig([.1 .1 1 .4],fontSize,1);

for r = 1:length(ROIs)
    f(1) = niceFig([.1 .1 .8 .8],fontSize,1);
    
    ROInum = cellNum(ROIs{r},info.ROIs);
    sFit = struct;
    for s = 1:length(subj)
        fits = subj(s).roi(ROInum).fits;
        
        for c = 1:length(fits)
            if c == 1 mult = .25; else mult = 1; end
            
            figure(f(1)); subplot(2,length(info.subjs)/2,s);
            
            x = [fits(c).vox.eccen]';
            y = [fits(c).vox.size]'/2;
            X = [x ones(length(fits(c).vox),1)];
            
            if ~isempty(X)
                [sFit(c).h(s,:),R2(s)] = fitl1line(X,y);
                h1 = plot(x,X*sFit(c).h(s,:)','Color',condColors(s,1)*mult); hold on;
                hold on; scatter(x,y,5,condColors(s,1)*mult,'filled'); hold on; title([ROIs{r} ' subj ' info.subjs{s}]);
                set(h1,'LineWidth',1); %set(extLine,'LineStyle',':');
                alpha(.2);
            end
        end
        figure(f(1));
        xl = xlim; xlim([0.25 6]); yl = ylim; if containsTxt(ROIs{r},'faces') ylim([0.25 6]); else ylim([0.25 6]);end
        set(gca,'TickDir','out');
        xlabel('Eccen (dva)','FontSize',fontSize); ylabel(['Size (1*SD/sqrt(N)) (dva)'],'FontSize',fontSize);
        axis square;
        
    end
    
    txt = ['fitl1line_' ROIs{r} '_' whichModel '_' whichStim '_r' num2str(minR2) '_' hemText(hems)];
    
    figure(f(2)); subplot(1,length(ROIs),r);
    for c = 1:2
        if c == 1 mult = .25; else mult = 1; end
        
        % full line fit
        x = [fitRange]';
        X = [x ones(length(x),1)];
        
        binThresh = 0;
        
        lines = sFit(c).h;
        [N,edges] = histcounts([roi(ROInum).fits(c).vox.eccen],fitRange);
        subFit = fitRange(N>binThresh)';
        fitOut = [subFit  ones(length(subFit),1) ] *lines';
        
        h1 = plot(fitRange,X*mean(lines)','Color',roiColors(ROIs{r})*c*.5); hold on;
        
        % construct errorbar
        
        plotErr = NaN(1,length(fitRange));
        plotErr(:,find(N>binThresh)) = se(fitOut');
        
        hold on; e = errorbar3(fitRange,[X*mean(lines)']',plotErr,'v',roiColors(ROIs{r})*c*.5); set(e,'FaceAlpha',.2);
        axis square;  xlim([0 6]);  ylim([0 6]); set(gca,'box','off','TickDir','out');
        title(ROIs{r});
    end
    
    %%% t test on line slope and intercept
    fprintf('---\n%s ttest:\n---\n',ROIs{r});
    
    
    tests = {'Slope' 'Intercept'};
    for t = 1:length(tests)
        [H,P,CI,STATS] = ttest(sFit(1).h(:,t),sFit(2).h(:,t));
        fprintf('%s: %s mean = %.3f (se = %.3f), %s mean = %.3f (se = %.3f)\n',tests{t},fits(1).cond, mean(sFit(1).h(:,t)),se(sFit(1).h(:,t)),...
            fits(2).cond, mean(sFit(2).h(:,t)),se(sFit(2).h(:,t)));
        if H addTxt = '***'; else addTxt = '';end
        fprintf('%st(%d) = %.3f, p = %.3f\n',addTxt,STATS.df, STATS.tstat,P);
    end
    
end

