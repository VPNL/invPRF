% coverage metrics (area, X, Y) for bilat visual regions

clear all; close all;

subjs = prfSubjs;
expt = 'fixPRF';


minR2 = 20;          % cutoff for vox selection
whichANOVA = 'face';
ROIs= [standardROIs(whichANOVA)];

whichStim = 'outline';
whichModel = 'kayCSS';
hems = {'rh' 'lh'};

tests = {'centY' 'centX' 'area'};

%%% how were the bootstraps generated?
boot.iters = 1000;
boot.vox = 0.8; % now implements this as a proportion of total voxels, not an absolute number
boot.method = 'binary';%'mean';%'max';%'binary';% % 'mean' or 'max'
boot.scaleSubjs = 0; % rescale each individual's coverage to [0 1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prfSet = ['prfSets/fixPRF_kayCSS_outline_' hemText(hems) '_r2-20.mat'];%pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems);
covPath = [pwd '/coverage/']; 

bootOpts = [boot.method '_iters' num2str(boot.iters) '_vox' num2str(boot.vox) '_scale' num2str(boot.scaleSubjs) '.mat'];
grpCov = [covPath fileName(prfSet) '_' bootOpts];
load(grpCov);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to run stats on?                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subjNum = cellNum(subjs,prfSubjs);

saveDir = [pwd '/stats/'];
checkDir(saveDir);
fid = fopen([saveDir 'ANOVA2_' hemText(hems) '_centroid_' whichANOVA '-' num2str(minR2) '_' whichModel '_' whichStim '.txt'],'w+');


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for r = 1:length(ROIs)
        
        ROInum = cellNum(ROIs{r},standardROIs);
        
        % centroid calc
            for c = 1:2
                for s = 1:length(subjNum)
                    [boot.roi(ROInum).cond(c).centX(s),cY] = FWHMcentroid(squeeze(boot.roi(ROInum).cond(c).cov(s,:,:)),.5);
                    boot.roi(ROInum).cond(c).centY(s) = boot.res+cY;
                end
 
            end
        end

        
   for t = 1:length(tests)
    anovaData= []; rmSubjs = [];
    fprintf(fid,'\n**************\n%s\n**************\n',strTogether(ROIs));
    fprintf(fid,'coverage file: %s\n',grpCov);     

    % aggregate data, check for missing values
     for r = 1:length(ROIs)
        ROInum = cellNum(ROIs{r},standardROIs);
        for s = 1:length(subjs)
            for c = 1:2
                eval(['sD = boot.roi(ROInum).cond(c).' tests{t} '(s);']);
                if ~isnan(sD)
                    anovaData = [anovaData; sD r c s];
                else rmSubjs(end+1) = s;
                end
            end
        end
    end
    
    fprintf('\n**************\n%s\n**************\n', strTogether(ROIs));
    
    % check for missing data and remove those subjects from the comparison
    [anovaData] = anova_rmSubjs(anovaData,rmSubjs);
    
    % run anova
    factNames = {'ROI' 'condition'};
    result = rm_anova2(anovaData(:,1),anovaData(:,end),anovaData(:,2),anovaData(:,3),factNames);
    %print output
    anova2_text(fid,result,tests{t});
     
    fprintf('\n*****\nt-tests by condition:\n****\n');
    %%% ROI-wise ttesting
    for r = 1:length(ROIs)
        rInd = anovaData(:,2) == r;
        rData = anovaData(anovaData(:,2) == r,:);
        for c = 1:2
            sData(:,c) = rData(rData(:,3)==c,1); end
        [H,P,CI,STATS] = ttest(sData(:,1),sData(:,2));
        diff = mean(sData(:,2)-sData(:,1));
        se = std(sData(:,2)-sData(:,1))/sqrt(length(sData));
        if strcmp(tests{t},'X') || strcmp(tests{t},'Y')
            diff = diff /roi(1).fits(1).ppd;
            se = se /roi(1).fits(1).ppd;end
        p(r) = P;
        stats{r} = STATS;
        if H sig = '***'; else sig = ''; end
        fprintf('%s [%s %s:] %s: t(%d)=%.5f, p=%.5f; diff = %.3f (SE=%.3f)\n',sig,hemText(hems),ROIs{r},tests{t},STATS.df,STATS.tstat,P,diff,se);
        fprintf(fid,'%s [%s %s:]  %s: t(%d)=%.5f, p=%.5f; diff = %.3f (SE=%.3f)\n',sig,hemText(hems),ROIs{r},tests{t},STATS.df,STATS.tstat,P,diff,se);
    end   
end

%playSound;