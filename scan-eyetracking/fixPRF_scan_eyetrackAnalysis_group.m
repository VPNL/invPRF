% analysis script for initial import and reading of eyelink files
% updated from prfRec2 to fixPRF5 scan

%%% ET NAME: always in the format subj/fixPRF_runNum - will be parsed later
clear all; close all;

subjs = prfSubjs;
saveFig = 0;
whichPre = '3'; % leave blank for no detrending, 2 for linear, 3 for 2nd order
if ~isempty(whichPre) preText = ['_pre' num2str(whichPre)]; else preText = ''; end


%s=2;
matDir = [pwd '/mats/'];

groupc = struct('name',{'blank' 'inverted' 'upright'},'samples',[],'subjs',[],'stdX',[],'stdY',[],'stdXY',[]);
for s = 1:length(subjs)
    %data from this subject
    
    subjFile = [matDir subjs{s} '/fixPRF_cdata_' subjs{s} preText '.mat'];% contains: 'subjc','ppd','centerPos'
    if exist(subjFile)
    fprintf('Loading %s...\n',subjFile);
    load(subjFile);
    
    for n = 1:length(subjc)
       groupc(n).name = subjc(n).name;
       groupc(n).samples = [groupc(n).samples;subjc(n).samples];
       groupc(n).subjs = [groupc(n).subjs subjs(s)];
       groupc(n).stdX = [groupc(n).stdX;subjc(n).stdX];
       groupc(n).stdY = [groupc(n).stdY;subjc(n).stdY];
       groupc(n).stdXY = [groupc(n).stdXY;subjc(n).stdXY];
    end
    else
        fprintf('*** Missing %s...\n',subjs{s});
    end 
end
               
    % main figure - collapsed across subjs
    niceFig([.1 .1 .8 .5]);
    for n = 1:length(groupc)
        subplot(1,4,n);
        groupc(n).convertedSamples = eyeInSpace_scan(subjc(n).samples,ppd,0.2,centerPos,5,1); % outputs a centered samples in dva
        title(groupc(n).name);
    end
    
    anovaData = []; colors = []; 
    for c = 1:length(groupc) % fact2
        for xy = 1:2 % fact1
            for s = 1:length(groupc(c).subjs)
            anovaData = [anovaData; groupc(c).stdXY(s,xy) c xy s];
            end
        end
        colors = [colors; condColors(c); lighter(condColors(c))];
    end
    
    %%% ANOVA
    factNames = {'X/Y' 'stim condition'};
    result = rm_anova2(anovaData(:,1),anovaData(:,end),anovaData(:,2),anovaData(:,3),factNames);
    %print output
    fid = fopen([pwd '/ANOVA2_N' num2str(s) preText '.txt'],'w+');
    anova2_text(fid,result,'std of eye movement');
    
   
    subplot(1,4,4);
     h= niceBars2(horzcat(groupc.stdXY),'mean',1,{'blankX' 'blankY' 'invX','invY','uprX','uprY'},colors);
    axis square;
    superTitle(['fixPRF N = ' num2str(s) preText]);
    
    if saveFig
        niceSave([pwd '/figures/group_results/'],['fixPRF_N' num2str(length(groupc(n).subjs)) preText ],{'png' 'svg'});
    end
%if onLaptop playSound; end
%end
checkDir([matDir 'group/']);
save([matDir 'group/groupData' preText '.mat'],'groupc','anovaData');