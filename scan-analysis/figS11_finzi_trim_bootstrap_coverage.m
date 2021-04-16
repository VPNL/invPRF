% code to trim & reanalyze finzi 2021 data to match poltoratski 2021
% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
prfSet = 'finzi2021_bi_pRFset_ve20.mat';
load(['prfSets/' prfSet]);

minR2 = 20; % default is 20 in this set - to further trim, edit this
ROIs= {'IOG_faces', 'pFus_faces', 'mFus_faces', 'pSTS_faces'};

computeCoverage = 0; % re-running this will generate a new bootstrap sample & will change downstream output!
contourPlot = 0; % cov plot or contour plot
plotTrims = 0;
plotIndivs =0;
plotCov = 1;
plotParStuff = 0;

boot.maxEccen = 5;   % in deg, max eccentricity. 50 in full set, 5 in limited set
boot.doTrim = 0;

% path to the manuscript's bootstrap samples
grpCov = ['coverage/finzi_eccen' num2str(boot.maxEccen)  '_bi_pRFset_ve20.mat'];
grpTrimmed = ['prfSets/trimmed_finzi_eccen' num2str(boot.maxEccen)  '_bi_pRFset_ve20.mat'];

plotCent = 1; % weighted centroid of each image
hem = 'bi';
boot.iters = 1000;
boot.vox = 0.8; % now implements this as a proportion of total voxels, not an absolute number
boot.method = 'binary';% 'max';%% 'mean' or 'max' or 'binary'
boot.szMult = 1; % if binary, factor by which size is multiplied at == radius. pre 5/4/20 this was 1, 2 matches kk and df
boot.scaleSubjs = 0; % rescale each individual's coverage to [0 1]
boot.subjs = info.subjs;
boot.hem = hem;
boot.minR2 = minR2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data & set up folders           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exptDir = [pwd '/finzi2021'];

if boot.szMult == 1 multText = ''; else multText = ['_mult' num2str(boot.szMult)]; end
bootOpts = ['analyzePRF_' boot.method '_minR2-' num2str(minR2)];
matPath = [exptDir '/'  bootOpts '/mats/eccen-' num2str(boot.maxEccen)]; checkDir(matPath);
covImPath = [exptDir '/'  bootOpts '/covPlots/eccen-' num2str(boot.maxEccen)]; checkDir(covImPath);
tic


if  ~exist(grpTrimmed) || computeCoverage
    allR = []; allC = [];
    
    %%% TRIMMING STEP
    
    for s = 1:length(info.subjs)
        for r =  1:length(ROIs)
            %%% get rid of condition spec
            subj(s).roi(r).vox = subj(s).roi(r).fits(1).vox;
            
            if ~isempty(subj(s).roi(r).vox) % if this subject has this ROI
                
                % par summary
                if plotTrims
                    niceFig([.1 .1 .5 .8]);
                    vox = subj(s).roi(r).vox;
                    subplot(2,2,1);
                    scatter([vox.eccen],[vox.size]); l= lsline; set(l,'lineWidth',1,'color','r');  axis square;
                    ylabel('size'); xlabel('eccen'); title([info.subjs{s} ' - no trim']);
                    subplot(2,2,2);
                    niceHist([vox.r2],condColors(1),1); title([info.subjs{s} ' r2 - no trim']); end
                %%% TRIM VOX TO MATCH INVPRF
                trim = [];
                if boot.doTrim
                for v = 1:length(subj(s).roi(r).vox)
                    vox = subj(s).roi(r).vox(v);
                    if vox.size > boot.maxEccen*2 || ...
                            vox.size < .1 || ...
                            vox.eccen > boot.maxEccen || ...
                            vox.r2 < minR2/100 || ...
                            isnan(vox.r2)
                        trim(end+1) = v;
                    end
                end
                else v = length(subj(s).roi(r).vox);
                end
                %%% try to find values that the algorithm has gotten stuck
                %%% on - this was an issue in original mrvista
                %%% implementation, but is resolved in analyzePRF
                %%% (finzi2020_bi_prfset) model fitting. analyzePRF #1!
                [reps,counts,inds]=findRepeats([subj(s).roi(r).vox.size],50); % cutoff = 50 reps/roi/subject in original prfset
                if ~isempty(reps)
                    for n = 1:length(reps)
                        fprintf('Size value %.3f has %i repetitions...\n',reps(n),counts(n)); end
                end
                allR = [allR reps]; allC = [allC counts];
                
                %%% trim those
                trim = [trim inds]; trim = unique(trim);
                
                subj(s).roi(r).vox(trim) = [];
                fprintf('Subj %s: trimmed %i of %i %s voxels...\n',info.subjs{s},length(trim),v,ROIs{r});
                if plotTrims
                    vox = subj(s).roi(r).vox;
                    subplot(2,2,3);
                    scatter([vox.eccen],[vox.size]); l= lsline; set(l,'lineWidth',1,'color','r');  axis square;
                    ylabel('size'); xlabel('eccen'); title([info.subjs{s} ' - post trim']);
                    subplot(2,2,4);
                    niceHist([vox.r2],condColors(1),1); title([info.subjs{s} ' r2 - post trim']);
                    superTitle([hem ' ' ROIs{r}],12,.05);s
                    niceSave([exptDir '/' bootOpts '/trimming/'],[info.subjs{s} '_' hem '-' ROIs{r} '_maxEccen' num2str(boot.maxEccen)]);
                    close(gcf);
                end
            else
                fprintf('** No voxels to bootstrap: %s %s%s .\n',info.subjs{s},hem,ROIs{r});
            end
        end
        
    end
    % to look at the repeated vals
    sizeRepeats = sortrows([allR;allC]',1,'descend');
    save(grpTrimmed,'subj','info');
    
else load(grpTrimmed);
end

if ~exist(grpCov) || computeCoverage
    for s = 1:length(info.subjs)
        % compute or load this subject's bootstrapped coverage
        
        for r =  1:length(ROIs)
            if ~isempty(subj(s).roi(r).vox)
                %%% MAIN BOOTSTRAPPING STEP
                sBoot(s).roi(r) = bootCoverage(subj(s).roi(r).vox,boot.method,boot.iters,boot.vox,boot.scaleSubjs,boot.szMult,boot.maxEccen*2);
                
                % centroid of FWHM
                [boot.roi(r).centX,boot.roi(r).centY] = FWHMcentroid(sBoot(s).roi(r).covIm);
                
                boot.roi(r).cov(s,:,:) = flipud(sBoot(s).roi(r).covIm);
                boot.roi(r).area(s) = sBoot(s).roi(r).areaDeg;
                boot.res = sBoot(s).roi(r).res; boot.ppd = sBoot(s).roi(r).ppd; % grab these, they are set within the bootstrap functo
                fprintf('Done with %s %s bootstrap step.\n',info.subjs{s},ROIs{r});
            else
                boot.roi(r).cov(s,:,:) = nan(111,111);
                fprintf('** No voxels to bootstrap: %s %s .\n',info.subjs{s},ROIs{r});
            end
            
        end
        
    end
    toc
    save(grpCov,'boot','sBoot');
    
else tic; load(grpCov); toc; end

if plotCov
    if plotIndivs
        for r = 1:length(ROIs)
            % plot individual coverage
            niceFig([.1 .1 1 .4]);
            for s = 1:length(boot.subjs)
                subplot(1,length(boot.subjs),s);
                if ~isnan(boot.roi(r).cov(s,1,1))
                    if s < length(boot.subjs) plotCovIm(squeeze(boot.roi(r).cov(s,:,:)),boot.res,boot.ppd,0,plotCent);
                    else plotCovIm(squeeze(boot.roi(r).cov(s,:,:)),boot.res,boot.ppd,1,plotCent);end
                    t = title({info.subjs{s}; ['FWHM: ' num2str(boot.roi(r).area(s)) 'dva']; ['Num Vox: ' num2str(length(subj(s).roi(r).vox))]}); set(t,'visible','on');
                else
                    axis square; set(gca,'visible','off');
                    t = title({info.subjs{s}; ['No Voxels']}); set(t,'visible','on');
                end
            end
            superTitle([hem '-' ROIs{r} ' ' fileName(prfSet)],14, .05);
            %niceSave(covImPath,[hem '-' ROIs{r}  '_'  fileName(grpCov) '.png']);
        end
        
    end
    
    meds = [];
    niceFig([.1 .1 .2*length(ROIs) .4]);
    for r = 1:length(ROIs)
        subplot(1,length(ROIs),r);
        
        % remove subjects who have insufficient (<10 voxels)
        for s = 1:length(boot.subjs)
            if length(subj(s).roi(r).vox) < 10
                boot.roi(r).cov(s,:,:) = nan(111,111);
                fprintf('Removing subject %s %s...\n',info.subjs{s},ROIs{r});
            end
        end
        
        % take overall coverage image
        meanIm = squeeze(nanmean(boot.roi(r).cov));
        
        if strcmp(boot.method,'max')
            meanIm(find(meanIm<.5)) = 0;
            %meanIm = imgaussfilt(meanIm,boot.ppd/10);
        end
        
        try
            if contourPlot plotCovContour(meanIm,boot.res,boot.ppd,1); else
                plotCovIm(meanIm,boot.res,boot.ppd,0,plotCent); end
            x = xlabel({'Mean FWHM:';sprintf('%.2f dva',mean(boot.roi(r).area));...
                sprintf('SE (N=%d) = %.2f',length(boot.subjs),se(boot.roi(r).area))});
            set(x,'visible','on');
        catch
            fprintf('Couldn''t plotCovIm %s %s...\n',hem,ROIs{r});
        end
        
        t = title(ROIs{r});
        set(t,'visible','on');
        
        superTitle(fileName(grpCov),14,.05);
        %niceSave(covImPath,[hem '-allROIs_eccen' num2str(boot.maxEccen) '_' fileName(grpCov)],[],[],{'png' 'svg'})
    end
    if onLaptop playSound; end
    
    
    for r = 1:length(ROIs)
        for s = 1:length(boot.subjs)
            meds(s,:) = nanmedian(reshape([subj(s).roi(r).vox.XYdeg],2,length(subj(s).roi(r).vox))');
            if ~isnan(boot.roi(r).cov(s,1,1))
                [centr(s,1),centr(s,2)] = FWHMcentroid(squeeze(boot.roi(r).cov(s,:,:)));
            else centr(s,:) = [nan nan]; end
        end
        centr(:,1);
        fprintf('%s -- Mean X: %.3f (SE=%.3f), Mean Y: %.3f (SE=%.3f)\n',ROIs{r},nanmean(meds(:,1)),se(meds(:,1)),nanmean(meds(:,2)),se(meds(:,2)));
        fprintf('%s -- Centroid X: %.3f (SE=%.3f), Centroid Y: %.3f (SE=%.3f)\n',ROIs{r},nanmean(centr(:,1)),se(centr(:,1)),nanmean(centr(:,2)),se(centr(:,2)));
        
    end
end

% plots other info from trimmed subj(s) file

% load('fixPRF_kayCSS_outline_bilat_r2-20.mat')

if plotParStuff
    miscPath = [exptDir '/'  bootOpts '/misc/eccen-' num2str(boot.maxEccen)];
    pos = {'X' 'Y'};
    whichPlot = 'summary';% 'hist'
    whichM = 'mean';
    
    switch whichPlot
        case 'hist'
            for t = 1:length(pos)
                for r = 1:length(ROIs)
                    % plot individual coverage
                    niceFig([.1 .1 1 .6]);sp = 1;
                    for s = 1:length(info.subjs)
                        subplot(2,ceil(length(info.subjs)/2),sp)
                        %XY= reshape([subj(s).roi(r).fits(2).vox.XYdeg],2,length(subj(s).roi(r).fits(2).vox))';
                        XY= reshape([subj(s).roi(r).vox.XYdeg],2,length(subj(s).roi(r).vox))';
                        
                        niceHist(XY(:,t),[condColors(1);condColors(2)],1)
                        title(subj(s).roi(r).fits(1).session)
                        xlabel(pos{t});
                        sp=sp+1;
                    end
                    superTitle([pos{t} ' ' hem '-' ROIs{r} ' ' fileName(prfSet)],14, .05);
                    %niceSave(miscPath,['indiv-hist_position-' pos{t} '_' hem '-' ROIs{r}  '_'  fileName(grpCov)]);
                end
                
            end
        case 'summary'
            niceFig([.1 .1 .2 .15*length(pos)]);
                dotColors = [];
                for t = 1:length(pos)
                    subplot(length(pos),1,t)
                    
                    for r = 1:length(ROIs)
                        for s = 1:length(subj)
                            XY= reshape([subj(s).roi(r).vox.XYdeg],2,length(subj(s).roi(r).vox))';
                            try
                                eval(['allPars(r,s) = nan' whichM '(XY(:,t));']);
                            catch allPars(r,s) = NaN; end % for missing values
                        end
                        dotColors = [roiColors(ROIs{r}); dotColors];
                    end
                    if strcmp(pos{t},'Y') == 1
                        allPars = allPars*-1; end % y flip
                    niceBoxPlusGrouped(allPars',ROIs,{['eccen-' num2str(boot.maxEccen)]},flipud(dotColors),[],1,0);
                    title(pos{t})
                end
            end
            superTitle([hem '-' ROIs{r} ' ' fileName(prfSet)],14, .05);
            %niceSave(miscPath,['summary-dots_position-' pos{t} '_' hem '-' ROIs{r}  '_'  fileName(grpCov)]);
                
    end
    playSound;