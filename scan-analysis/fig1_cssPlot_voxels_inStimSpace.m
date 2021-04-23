% loads & plots some basics re: invPRF fits (css model)
% SP 10/25/17
% _indiv code allows us to plot individual ROIs in multiple plots
% 10/8/18 edit: made more flexible across expriment with different numbers
% of fit conditions

clear all; close all;

% ROIs and subjects
subjs = prfSubjs;
expt = 'fixPRF';
ROIs= standardROIs([1 7]);%'face');%{'V1'};%['hV4' standardROIs('face')];% %{'PL' 'ML'};%

% what to plot
plotFit = 1; % plot the actual voxel fit as well as the beta rings
sortR2 = 1; % if 0, grab random voxels
sampleVox = 50; % how many randomly selected voxels are we plotting?

% plotting choices
saveFig = 1;
im.gridSpaceDeg = 1.5;
im.faceSizeDeg = 3.2;
im.ppd = 100; % these images need to be fairly high-res


% fit file
minR2 =20;             % cutoff for vox selection
whichStim = 'outline'; %
whichModel = 'kayCSS'; %
hems = {'rh' 'lh'};

plotOrder = [2 1];
plotX = [1:25; 28:52];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 12; titleSize = 14;

% now we load in the data from both hemispheres, and threshold across
load(['prfSets/fixPRF_kayCSS_outline_' hemText(hems) '_r2-20.mat']);
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

for r = 1:length(ROIs)
    if length(subjNum) == 1 fits = subj(subjNum).roi(ROInum(r)).fits;
    else fits =roi(ROInum(r)).fits; end
    
    [plotVox] = pickPlotVox(fits,sortR2,sampleVox);
    
    for vv = 1:length(plotVox)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        niceFig([0 .2 1 .5],fontSize);
        numPlots = [1 5];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        v = plotVox(vv);
        cl = [min(min([fits(1).vox(v).betas fits(2).vox(v).betas]))...
            max(max([fits(1).vox(v).betas fits(2).vox(v).betas]))];
        subplot(numPlots(1),numPlots(2),1:2);
        
        for c = 1:length(fits)
            hbar(c) = niceBars_colValue(plotX(c,:),[fits(plotOrder(c)).vox(v).betas],[fits(plotOrder(c)).vox(v).sems],cl(1),cl(2),colormap('parula'),.6);
            plot(plotX(c,:),fits(plotOrder(c)).vox(v).modelfit,'k-','LineWidth',.5);
        end
        
        xlabel(strTogether({fits(plotOrder).cond},10)); xlim([0 plotX(2,end)]); set(gca,'xticklabel',[]);
        ylabel('Beta Estim.'); pbaspect([2,.5,1]);%axis square;
        
        for tt = 1:length(fits)
            parText{tt} = ['R^{2}=' num2str(round(fits(tt).vox(v).r2)) ' ' fits(tt).cond ];
        end
        t = title(parText);
        
        set(gca,'box','off','FontName','Arial','FontWeight','normal','TickDir','out');
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % beta value plots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for c = 1:length(fits)
            
            subplot(numPlots(1),numPlots(2),2+c)
            betas = reshape(fits(plotOrder(c)).vox(v).betas,sqrt(length(fits(plotOrder(c)).vox(v).betas)),sqrt(length(fits(plotOrder(c)).vox(v).betas)))';
            
            betasInStim(betas,im.ppd,im.gridSpaceDeg,im.faceSizeDeg,cl,1); axis square;
            
            t = title([fits(plotOrder(c)).cond ' Betas']);
            set(t,'visible','on');
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% plot prf estimates
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(numPlots(1),numPlots(2),numPlots(2));
        for c = 1:length(fits)
            % PRF size is defined as S/sqrt(N).
            h(c) = plotCoverage(fits(c).vox(v),condColors(c),[],fits(1).ppd,fits(1).res,1,0);
            if c == length(fits) l=legend([h(1:length(fits))],{fits.cond});
                set(l,'box','off','location','Best');
                for tt = 1:length(fits)
                    parText{tt} = ['Params' num2str(tt) ': [' num2str(round(fits(tt).vox(v).params)) ']'];
                end
                t = title(parText);
                set(t,'visible','on');
            end
        end
        
        titleText = [ROIs{r} ' Voxel #' num2str(v) ', ' whichStim ' stim, ' whichModel ' model'];
        
        titleText = [titleText ', Session = ' fileName(fits(1).vox(v).stim,8)];
        superTitle(titleText,12, .05)
        if saveFig
            txt = [whichModel  '_' whichStim '-r' num2str(vv) '-vox' num2str(v)];
            niceSave([pwd '/voxPlots/' hemText(hems) '_' ROIs{r} '/'],txt,[],[],{'png'}); % just save pngs, since these can be generated pretty quickly
            if sampleVox > 10; close(gcf); end
        end
    end % vox
end
if onLaptop playSound; end