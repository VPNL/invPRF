% primary analysis script for behavioral experiment
clear all; close all;

% establish directories
exptDir = [dirOf(pwd) 'behavior-experiment/']; 
dataDir = [exptDir 'data/'];

% load log file
load([exptDir 'prfRec2_log.mat'])

an.subjs = {session.name};
an.expt = 'prfRec2';
an.exclude = {'MH' 'MG' 'JJ'};
an.which = 'scan';  % scan or all, for now (can also add lab or sona or other distinctions below)
an.plotSubjs = 0;   % binary, plot individual subject performance
saveFig = 0;
doRemove = 0;   % re-run the trial removal process

%%% choose subjects accordingly
switch an.which
    case 'scan'
        an.subjs = intersect(prfSubjs,an.subjs);
        an.subjs = setdiff(an.subjs,an.exclude);
    case 'all'
        an.subjs = setdiff(an.subjs,an.exclude);
end

if an.plotSubjs [sd1,sd2] = subplotDims(length(an.subjs)); f(1) = niceFig([.1 .1 1 1]); f(2) = niceFig([.1 .1 1 1]);end



% aggregate
for s = 1:length(an.subjs)
    if isfile([dataDir 'prfRec2_' an.subjs{s} '_withQA.mat']) && ~doRemove
        load([dataDir 'prfRec2_' an.subjs{s} '_withQA.mat']);
    else
    df = [dataDir 'prfRec2_' an.subjs{s} '.mat'];
    try
        load(df); fprintf('Loading subject file: %s...\n',df);
    catch
        fprintf('Error loading subject file: %s.\n',df);
    end
    
    try
    % load QA xls data
    qa = xlsread([dirOf(exptDir,1),'analysis/eyetrackQA.xlsx'],an.subjs{s},'B2:S21');
    catch
        fprintf('Error loading subject QA file: %s\n',an.subjs{s});
    end
    removed = [];
    
    % remove trials on which subject broke fixation
    for b = 1:length(block)
        rem = find(qa(:,b)==0);
        removed = [removed block(b).trial(rem).cond];
                block(b).trial(rem) = []; end
            
    fprintf('*** Subj: %s ***\n',an.subjs{s});
    for c = 1:length(condition)
       an.removed(s,c) = length(find(removed==c))/params.condTrials;
       fprintf('%0.2f%% trials removed: %s.\n',100*an.removed(s,c),condition(c).name);
    end
  
 
% performance
perf = struct('condName',{condition.name},'hits',[],'CRs',[],'hitRate',[],'CRrate',[]);

for c = 1:length(perf)
    for b = 1:length(block)
        ind = find([block(b).trial.cond]==c); % find trials corresponding to this condition
        same = ind(find([block(b).trial(ind).targ]==1)); % of those, find same trials
        diffr = ind(find([block(b).trial(ind).targ]==0)); % and different trials
        perf(c).hits = [perf(c).hits block(b).trial(same).correct];
        perf(c).CRs = [perf(c).CRs [block(b).trial(diffr).correct]];
    end
    perf(c).CRrate = nanmean([perf(c).CRs]);
    perf(c).hitRate = nanmean([perf(c).hits]);
    perf(c).percCorrect = (perf(c).hitRate + perf(c).CRrate)/2 * 100;

    
    %%% do we need to apply correction to this dprime rate?
    % correct for 0 false alarms by setting it to 0.5 false alarms, and 1
    % hit rate by setting it to 0.5 errors 
    
    
    txt = 'uncorrected';
    if perf(c).CRrate == 1 FAfix = .25/params.condTrials; txt = '*corrected'; else FAfix = 1-perf(c).CRrate; end
    if perf(c).hitRate == 1 HRfix = 1-.25/params.condTrials; txt = '*corrected'; else HRfix = perf(c).hitRate; end
    
    perf(c).dprime = norminv(HRfix)-norminv(FAfix);
    
    fprintf('%s Condition: %0.2f%% correct responses, %s d'' = %0.2f.\n',perf(c).condName,perf(c).percCorrect,txt,perf(c).dprime);
end
    % save a version of this performance
    save([dataDir 'prfRec2_' an.subjs{s} '_withQA.mat'],'block','condition','el','expt','params','perf');
    end
    
    % aggregate over subjects
    if s==1 grp = struct('cond',{perf.condName},'subjs',[],'subjPC',[],'subjD',[]); end % initalize the group struct
    for c = 1:length(grp)
        grp(c).subjs = [grp(c).subjs an.subjs(s)];
        grp(c).subjPC = [grp(c).subjPC perf(c).percCorrect];
        grp(c).subjD = [grp(c).subjD perf(c).dprime];
    end
    
    if an.plotSubjs
        figure(f(1));
        subplot(sd1,sd2,s);
        prfRec_plot([perf.percCorrect],'Percent Correct',{'Low' 'Center' 'Up'},[4 5 6 1 2 3]); title(an.subjs{s});
        figure(f(2));
        subplot(sd1,sd2,s);
        prfRec_plot([perf.dprime],'d Prime',{'Low' 'Center' 'Up'},[4 5 6 1 2 3]);title(an.subjs{s}); end
end

for c = 1:length(grp)
    grp(c).PCmean = nanmean(grp(c).subjPC);
    grp(c).PCse = nanstd(grp(c).subjPC)/sqrt(length(grp(c).subjPC));
    grp(c).Dmean = nanmean(grp(c).subjD);
    grp(c).Dse = nanstd(grp(c).subjD)/sqrt(length(grp(c).subjD));
end

save([nStr(an.subjs) '-' an.expt '.mat'],'grp','an');
    
grpNames = {'LowerLeft' 'Center' 'UpperRight'}; barNames = {'Upright' 'Inverted'};
order = [4 1; 5 2; 6,3];

metrics = {'dPrime' 'Percent Correct'}; niceFig([.1 .1 .6 .6]);
for n = 1:length(metrics)
    if strcmp(metrics{n},'dPrime') data = [grp.Dmean]; err = [grp.Dse];
    else data = [grp.PCmean]; err = [grp.PCse]; end
    subplot(1,2,n)
    niceGroupedBars(data(order),err(order),grpNames,barNames);
    title([metrics{n} ',  N = ' num2str(length(an.subjs)) ', Subjs: ' strTogether(an.subjs)]);
    ylabel(metrics{n});
end
if saveFig niceSave([pwd '/figures/'],[an.expt '_' nStr(an.subjs) '_groupPerf_withQA']);
    if an.plotSubjs
figure(f(1)); niceSave([pwd '/figures/'],[an.expt '_' nStr(an.subjs) '_indiv_' metrics{1} '_withQA']);
figure(f(2)); niceSave([pwd '/figures/'],[an.expt '_' nStr(an.subjs) '_indiv_' metrics{2} '_withQA']);
    end
end
    
%%%% ANOVA TIME!
  
    % Parameters:
    %    Y          dependent variable (numeric) in a column vector
    %    S          grouping variable for SUBJECT
    %    F1         grouping variable for factor #1
    %    F2         grouping variable for factor #2
    %    FACTNAMES  a cell array w/ two char arrays: {'factor1', 'factor2'}
    %
    %    Y should be a 1-d column vector with all of your data (numeric).
    %    The grouping variables should also be 1-d numeric, each with same
    %    length as Y. Each entry in each of the grouping vectors indicates the
    %    level # (or subject #) of the corresponding entry in Y.
    %    stats = rm_anova2(Y,S,F1,F2,FACTNAMES)
    
    sn = length(grp(1).subjs); cn = length(grp);
    anovaData = horzcat(grp.subjPC)';
  
    
    subjInd= repmat([1:sn],1,cn)';
    inverted = kron([1:2],ones(1,sn*3))';
    location = repmat(kron([1:3],ones(1,sn)),1,2)';
    [anovaData,subjInd,inverted,location];
    
results = rm_anova2(anovaData,subjInd,inverted,location,{'Up/Inv','Location'})


FIE = anovaData(find(inverted==2)) - anovaData(find(inverted==1));
loc = location(find(inverted==1)); % just grabbing our location indexing for 1/2 the data
sI = subjInd(find(inverted==1));


for l = 1:3
[H,P,CI,stats] =ttest(FIE(find(loc==l)));%;,FIE(find(fact==3)));
fprintf('T-test on FIE at location %s: t(%d) = %.4f, p = %.4f\n',grpNames{l},stats.df,stats.tstat,P);
end
