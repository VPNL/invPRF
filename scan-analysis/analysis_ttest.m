% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs = prfSubjs;
expt = 'fixPRF';
tests = {'Ydeg' 'Xdeg' 'eccen' 'gain' 'size' 'r2'}; % can be parname (Y,X,sd,gain,exp,shift) or pRF.read value (r2,size,eccen,gain)
whichM = 'mean'; % mean or median


r2cutoff = 'r2-20';         % cutoff for vox selection
ROIs= ['V1' standardROIs('face')];%('face-');%


whichStim = 'outline';
whichModel = 'kayCSS';
hems = {'rh' 'lh'};
txtName = [hemText(hems) '_' r2cutoff];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['prfSets/fixPRF_kayCSS_outline_' hemText(hems) '_r2-20.mat']);
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);
checkDir([pwd '/stats']);
fid = fopen([pwd '/stats/ttests_' whichM '_' txtName '_' whichModel '_' whichStim '.txt'],'w+');


comps = nchoosek(1:length(roi(1).fits),2);
fprintf('\n%s %s\n\n**************\n',whichM, r2cutoff);
fprintf(fid,'\n%s %s\n\n**************\n',whichM,r2cutoff);
for t = 1:length(tests)
    test = tests{t};
    fprintf('-----\n%s:\n-----\n',test);
    fprintf(fid,'-----\n%s:\n-----\n',test);
    for cc = 1:size(comps,1)
        baseCond = comps(cc,1);
        c = comps(cc,2);
        comp(c).descr = [roi(1).fits(baseCond).cond ' vs ' roi(1).fits(c).cond];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  ttest the parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        parNames = roi(1).fits(1).parNames;
        testNum = cellNum(test,parNames);
        
        
        for r = 1:length(ROIs)
            sData = [];
            for s = 1:length(subjs)
                % grab basecond
                vox = subj(subjNum(s)).roi(ROInum(r)).fits(baseCond).vox;
                % grab comp cond
                vox2 = subj(subjNum(s)).roi(ROInum(r)).fits(c).vox;
                if ~isempty(vox)
                    if ~isempty(testNum)
                        pars = vertcat(vox.params);
                        sData(s,1) = nanmean(pars(:,testNum));
                        pars2 = vertcat(vox2.params);
                        sData(s,2) = nanmean(pars2(:,testNum));
                    else
                        eval(['sData(s,1) = nan' whichM '([vox.' test ']);']);
                        eval(['sData(s,2) = nan' whichM '([vox2.' test ']);']);
                    end
                else sData(s,1) = NaN; sData(s,2) = NaN;
                end
            end
            
            %%% accounting for missing data
            try
                sData = rmmissing(sData);
            catch warning('Can''t run rmmissing() on this version of Matlab!'); end
            
            comp(c).groupData{r} = sData;
            [H,P,CI,STATS] = ttest(sData(:,1),sData(:,2));
            comp(c).diff = mean(sData(:,2)-sData(:,1));
            comp(c).se = std(sData(:,2)-sData(:,1))/sqrt(length(sData));
            if strcmp(test,'X') || strcmp(test,'Y')
                comp(c).diff = comp(c).diff /roi(1).fits(1).ppd;
                comp(c).se = comp(c).se /roi(1).fits(1).ppd;end
            comp(c).p(r) = P;
            comp(c).stats{r} = STATS;
            if H sig = '***'; else sig = ''; end
            fprintf('%s [%s %s:] %s param %s, %s: t(%d)=%.2f, p=%.3f; diff = %.3f (SE=%.3f)\n',sig,hemText(hems),ROIs{r},test,comp(c).descr,whichM,STATS.df,STATS.tstat,P,comp(c).diff,comp(c).se);
            fprintf(fid,'%s [%s %s:] %s param %s, %s: t(%d)=%.2f, p=%.3f; diff = %.3f (SE=%.3f)\n',sig,hemText(hems),ROIs{r},test,comp(c).descr,whichM,STATS.df,STATS.tstat,P,comp(c).diff,comp(c).se);
            
        end
    end
end

if onLaptop playSound; end
