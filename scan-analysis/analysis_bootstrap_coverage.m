% creates and saves a bootstrap-ed coverage map for every individual
% subject, and then averages them together for the total map.
% hybrid of mrvista coverage generation (bootstrapping) and SP coverage
% generation (no spatial smoothing, other choices...)
% note - re-running this bootstrap may alter the numerical output of
% downstream analyses! the iteration used in the paper is found in
% coverage/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

subjs = prfSubjs;
expt = 'fixPRF';

minR2 = 'r2-20';        % cutoff for vox selection
ROIs= standardROIs;
plotIndivs = 1;
computeCoverage = 1;
contour = 0;    % plot image or contour
   

whichStim = 'photo';
whichModel = 'kayCSS';
plotCent = 1; % weighted centroid of each image
hems = {'lh' 'rh'};

boot.iters = 1000;
boot.vox = 0.8; % now implements this as a proportion of total voxels, not an absolute number
boot.method = 'binary';% 'max';%% 'mean' or 'max' or 'binary'
boot.szMult = 2; % if binary, factor by which size is multiplied at == diameter. pre 5/4/20 this was 1 (becuase 2xsize was set in readPRFs.m), 2 matches kk and df and is appropriate for all summer-2020+-generated pRFsets
boot.scaleSubjs = 0; % rescale each individual's coverage to [0 1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prfSet = ['prfSets/fixPRF_kayCSS_outline_' hemText(hems) '_r2-20.mat'];
load(prfSet);

ROInum = cellNum(ROIs,info.ROIs)';
subjNum = cellNum(subjs,info.subjs);

exptDir = pwd;
covPath = [exptDir '/coverage/']; checkDir(covPath);

if boot.szMult == 2 multText = ''; else multText = ['_mult' num2str(boot.szMult)]; end
bootOpts = [boot.method '_iters' num2str(boot.iters) '_vox' num2str(boot.vox) '_scale' num2str(boot.scaleSubjs) multText '.mat'];


tic
grpCov = [covPath fileName(prfSet) '_' bootOpts];
if ~exist(grpCov) || computeCoverage
    
    for s = 1:length(info.subjs)
        sCov = [covPath info.subjs{s} '_' fileName(prfSet) '_' bootOpts];
        % compute or load this subject's bootstrapped coverage
        if ~exist(sCov) || computeCoverage
            for r = 1:length(info.ROIs)
                for c = 1:length(subj(1).roi(1).fits)
                    if ~isempty(subj(s).roi(r).fits(c).vox) % if this subject has this ROI
                        %%% MAIN BOOTSTRAPPING STEP
                        sBoot.roi(r).cond(c) = bootCoverage(subj(s).roi(r).fits(c).vox,boot.method,boot.iters,boot.vox,boot.scaleSubjs,boot.szMult);
                        
                        % centroid of FWHM
                        [boot.roi(r).cond(c).centX,boot.roi(r).cond(c).centY] = FWHMcentroid(sBoot.roi(r).cond(c).covIm);
                        
                        boot.roi(r).cond(c).cov(s,:,:) = sBoot.roi(r).cond(c).covIm;
                        boot.roi(r).cond(c).area(s) = sBoot.roi(r).cond(c).areaDeg;
                        boot.res = sBoot.roi(r).cond(c).res; boot.ppd = sBoot.roi(r).cond(c).ppd; % grab these, they are set within the bootstrap functo
                    end
                end
            end
            
        else load(sCov); end
    end
    toc
    save(grpCov,'boot');
else tic; load(grpCov); toc; end

if plotIndivs
    for r = ROInum
        % plot individual coverage
        niceFig([.1 .1 1 .4]); sp=1;
        for c = 1:length(boot.roi(1).cond)
            for ss = 1:length(subjNum)
                s = subjNum(ss);
                subplot(length(boot.roi(1).cond),length(subjNum),sp); sp=sp+1;
                if ss < length(subjNum) plotCovIm(squeeze(boot.roi(r).cond(c).cov(s,:,:)),boot.res,boot.ppd,0,plotCent);
                else plotCovIm(squeeze(boot.roi(r).cond(c).cov(s,:,:)),boot.res,boot.ppd,1,plotCent);end
                if ss == 1 y = ylabel(roi(1).fits(c).cond,'FontSize',20); set(y,'visible','on');  end
                t = title({info.subjs{s}; ['FWHM: ' num2str(boot.roi(r).cond(c).area(s)) 'dva']}); set(t,'visible','on');
            end
        end
        superTitle([ROIs{find(ROInum==r)} ' ' fileName(prfSet)],14, .05);
        % niceSave([raid 'invPRF/figures/' expt '/coverage/bootstrap/indivs/'],[ROIs{r} '_' fileName(grpCov)]);
    end
    
end


for r = ROInum
    niceFig([.02*r .02*r .8 .5]);
    % plot mean coverage

    for c = 1:length(boot.roi(1).cond)
        subplot(1,length(boot.roi(1).cond),c);
        % take overall coverage image
         meanIm = squeeze(nanmean(boot.roi(r).cond(c).cov));
         
        if strcmp(boot.method,'max')
            meanIm(find(meanIm<.5)) = 0;
                %meanIm = imgaussfilt(meanIm,boot.ppd/10);
        end

        if contour plotCovContour(meanIm,boot.res,boot.ppd,1); else
            plotCovIm(meanIm,boot.res,boot.ppd,1,plotCent); end
        x = xlabel({'Mean FWHM:';sprintf('%.2f dva',mean(boot.roi(r).cond(c).area));...
            sprintf('SE (N=%d) = %.2f',length(subjNum),se(boot.roi(r).cond(c).area))});
        set(x,'visible','on');
        y= ylabel(roi(1).fits(c).cond,'FontSize',20); set(y,'visible','on');
        t = title({[ROIs{find(ROInum==r)}];sprintf('t(%d) = %.2f',STATS.df,STATS.tstat);sprintf('p = %.3f',P)});
        set(t,'visible','on');if H set(t,'Color',[0 .5 0]); end
    end
    superTitle(fileName(grpCov),14,.05);
    % niceSave([raid 'invPRF/figures/' expt '/coverage/bootstrap/'],[ROIs{r} '_' fileName(grpCov)],[],[],{'svg' 'png'});
end
if onLaptop playSound; end
