% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs = prfSubjs;
expt = 'fixPRF';

saveFig =1;

minR2 = 'r2-20';          % cutoff for vox selection
ROIs= ['V1' standardROIs('face')];%standardROIs;%

whichStim = 'outline';%
whichModel = 'kayCSS';%'
whichPlot = 'split'; %'line'; %'shaded'% or %'shaded' or 'bars'
yOrEccen = 'y';

hems = {'rh' 'lh'};
binWidth = .25; % histogram bins

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;


% now we load in the data from both hemispheres, and threshold across
load(['prfSets/fixPRF_kayCSS_outline_' hemText(hems) '_r2-20.mat']);
niceFig([.1 .1 .4 .3]);

for r = 1:length(ROIs)
    ROInum = cellNum(ROIs{r},info.ROIs);
    for c = 1:2
        
        counts = [];
        for s = 1:length(subjs)
            switch yOrEccen
                case 'eccen'
                    edges = [.5:binWidth:5.25]; yl = [0 .2]; xl = [0 5];
            counts(s,:) = histcounts([subj(s).roi(ROInum).fits(c).vox.eccen],edges,'Normalization','probability');
                case 'y'
                    edges = [-5:binWidth:5];yl = [0 .2]; xl = [-5 5];
           counts(s,:) = histcounts([subj(s).roi(ROInum).fits(c).vox.Ydeg],edges,'Normalization','probability');
            end
        end
        allCounts{r,c} = counts;
        switch whichPlot
            case 'shaded'
                subplot(1,length(ROIs),r)
                hold on; plot([edges(1:end-1)+binWidth/2],mean(counts),'Color',roiColors(ROIs{r})*(c*.5));
                hold on; e(c) = errorbar3([edges(1:end-1)+binWidth/2],mean(counts),se(counts),'v',roiColors(ROIs{r})*(c*.5)); set(e(c),'FaceAlpha',.5);
            case 'bars'
                %subplot(2,length(ROIs),length(ROIs)*(c-1)+r);
                subplot(1,length(ROIs),r)
                if c == 1; mult = .7; else mult = .3; end
                niceBars_colValue([edges(1:end-1)+(mult*binWidth)],mean(counts),se(counts),yl(1),yl(2),colormap('parula')); hold on;
            case 'split'
                subplot(2,length(ROIs),r+(length(ROIs)*(c-1)))
                if c == 1; mult = .7; else mult = .3; end
                niceBars_colValue([edges(1:end-1)],mean(counts),se(counts),yl(1),yl(2),colormap('parula'),binWidth+.25); hold on;

                xlim(xl); xticks([xl(1):1:xl(2)]); xlabel('bin centers');%xticklabels(edges(1:end-1));
                set(gca,'TickDir','out'); ylim(yl); set(gca, 'box','off'); pbaspect([2 1 1]); %axis square;
                if c == 1 title([ROIs{r}]); end

            case 'line'
                subplot(1,length(ROIs),r)
                hold on; errorbar([edges(1:end-1)],mean(counts),se(counts),'Color',roiColors(ROIs{r})*(c*.5));
                
        end
        if ~strcmp(whichPlot,'split')
            xticks([xl(1):1:xl(2)]); xlabel('bin centers');%xticklabels(edges(1:end-1));
            set(gca,'TickDir','out'); ylim(yl); set(gca, 'box','off');
            title([ROIs{r}]);
            xlabel(yOrEccen); ylabel('Proportion of pRF centers'); %legend(e,{roi(1).fits(1).cond})
            axis square;end
    end
    
end


titleText = [expt ' (' hemText(hems) '), voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
superTitle(titleText,titleSize,.025);

% if saveFig == 1
%     txt = [yOrEccen '_acrossROIs_mockup' whichPlot '_' hemText(hems)];
%     niceSave([dirOf(pwd) 'figures/' expt '/centersByEccen/'],txt,[],[],{'svg' 'png'}); % just save pngs, since these can be generated pretty quickly
% end

if onLaptop playSound; end