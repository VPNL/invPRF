% compares behavioral performace across locations to pRF coverage for
% selected ROIs

% SP 8/8/18
clear all; close all;
expt = 'fixPRF';
saveFig = 0;

subjs = prfSubjs;
pooled = 0; % plot not as individual  locations but as center/mid/inner rings - rough eccen pooling
rings = {[[1:5:21] [5:5:25] [2:4] [22:24]];...
    [[7:9] [17:19] [12,14]];...
    [13]}; % this isn't perfect in .pos values (since those are not indexed in native matlab convention,


grpPerf = [];

for s = 1:length(prfSubjs)

subj = subjs{s};
% load data
[session, numRuns] = vpnlSessions(expt,subj); % OPTIONAL: SESSNUM, TASK

% but since we're just making rings, it works

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  create supertitle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

titleText = [expt ' Behav Performance, Subj:' subj ', Session ' session ' (' num2str(numRuns) ' Runs)'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1:
% 1) letter task performance, inverted
% 2) letter task performance, upright
% 3) letter task performance, blank trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fontSize = 14; titleSize = 18;
niceFig([.1 .1 .8 .5],fontSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% task 1
switch expt% load subject's performance
    case 'fixPRF'
        [stim,perfPlot] = fixPRF_getBehavior(session,numRuns,[dirOf(pwd) expt],1);
    case 'compPRF'
        [stim,perfPlot] = compPRF_getBehavior(session,numRuns,[dirOf(pwd) expt],1);
    case 'invPRF3'
        [stim,perfPlot] = invPRF3_getBehavior(session,numRuns,[dirOf(pwd) expt],1);
end
numPlots = [1 size(perfPlot,2)];
for c = 1:length(perfPlot)
    subplot(numPlots(1),numPlots(2),c)
    if pooled == 1
        for b = 1:length(rings)
            perfPlot(c).mat(rings{b}) = nanmean(perfPlot(c).mat(rings{b}));
        end
    end
    plotInSpace(perfPlot(c).mat,'Behavior',[perfPlot(c).name ': ' num2str(nanmean(perfPlot(c).vect)) ' hit rate'],1,[0 1]);
    
    grpPerf(s,c) = nanmean(perfPlot(c).vect);
end

if pooled == 1
    titleText = ['POOLED by Rough Eccen, ' titleText]; end
superTitle(titleText,titleSize-4,.01)

if saveFig == 1
    if pooled == 1
        txt = ['pooledBehav_' session];
        niceSave([dirOf(pwd) 'figures/' expt '/scanBehav/'],txt); % just save pngs, since these can be generated pretty quickly
    else
        txt = ['behav_' session];
        niceSave([dirOf(pwd) 'figures/' expt '/scanBehav/'],txt); % just save pngs, since these can be generated pretty quickly
    end
end

end

p = mean(grpPerf,2);
fprintf('Overall performance: mean = %.2f (SE = %.2f), range = %.2f - %.2f\n',mean(p),se(p),min(p),max(p));

% difference between upright and inverted faces?

 [~,p,~,stats] = ttest(grpPerf(:,1),grpPerf(:,2))
 fprintf('Performance difference between upright and inverted faces: %.3f (se = %.3f) Upright %.3f (se = %.3f) Inverted, t(%d) = %.3f, p = %.3f)\n',mean(grpPerf(:,2)),se(grpPerf(:,2)),mean(grpPerf(:,1)),se(grpPerf(:,1)),stats.df,stats.tstat,p);
