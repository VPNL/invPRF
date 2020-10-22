% analyze prfRec data
clear all; close all;

% establish directories
if ~isdir('/Volumes/projects') exptDir = '/share/kalanit/biac2/kgs/projects/behavFIE/prfRec/'; else
    exptDir = '/Volumes/projects/behavFIE/prfRec/'; end
dataDir = [exptDir 'data/'];

% load log file
load([exptDir 'prfRec2_log.mat'])

an.subjs = {session.name};
an.expt = 'prfRec2';
an.exclude = {'MH' 'MG' 'JJ'};% {''};%
an.which = 'scan';  % scan or all, for now (can also add lab or sona or other distinctions below)
an.plotSubjs = 0;   % binary, plot individual subject performance
saveFig = 0;
doRemove = 0;   % re-run the trail removal process

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
end


order = [4 1 5 2 6 3]; labels = {'Upr-LL' 'Inv-LL' 'Upr-Cen' 'Inv-Cen' 'Upr-UR' 'Inv-UR'};%
colors = repmat([condColors(1);condColors(1).*.25],3,1);

metrics = {'dPrime'}; niceFig([.1 .1 .6 .6]);
for n = 1:length(metrics)
    if strcmp(metrics{n},'dPrime') data = vertcat(grp.subjD)';
    else data = vertcat(grp.subjPC)'; end
    data = data(:,order)
    [hbar] = niceBars2(data,'mean',3,[],colors)
    title([metrics{n} ',  N = ' num2str(length(an.subjs)) ', Subjs: ' strTogether(an.subjs)]);
    ylabel(metrics{n});
end
if saveFig niceSave([pwd '/figures/'],[an.expt '_' nStr(an.subjs) '_scatterSubjs']);
    
end
