% analysis script for in-scanner eyetracking analysis on individual
% subject. works with .mat files generated via
% fxPRF_scan_eyetrackPreprocessing_*. to save space, original asc and edf
% files are not included in this repo but can be provided on request.

%%% ET NAME: always in the format subj/fixPRF_runNum - will be parsed later
clear all; close all;

subjs = {'SP'};%prfSubjs;
whichPre = '3'; % leave blank for no detrending, 2 for linear, 3 for 2nd order
if ~isempty(whichPre) preText = ['_pre' num2str(whichPre)]; else preText = ''; end

runFigs = 1;
saveFigs = 0;

for s = 1:length(subjs)
    
    subjc = struct('name',{'blank' 'inverted' 'upright'},'samples',[],'runInd',[]);
    
    for r = 1:10
        % input
        etName = sprintf('%s/fixPRF_%i',subjs{s},r);
        
        %dataFile = [exptDir 'data/' dataName  '.mat'];
        matDir = [pwd '/mats/'];
        etMat =  [matDir etName '_preprocessed' whichPre '.mat'];
        scanMat = [dirOf(pwd) 'scan-experiment/run-output/' subjs{s} '/fixPRF_' num2str(r) '.mat'];
        figDir =  [pwd '/figures/'];
        
        if exist(etMat)
            fprintf('Starting %s...\n',etMat);
            load(scanMat); load(etMat);
            
            
            ppd = scan.ppd;
            fixRad = params.fixRadDeg;
            %eyeInit.screen.width = scan.screenWidth;
            
            % calculate center position as median across conditions
            centerPos = nanmedian(samples(:,2:3));
            
            
            if runFigs
                niceFig([.1 .1 .9 .6]);
                subplot(1,4,1);
                
                % all trial over time
                eyeInSpace_scan(samples,ppd,fixRad,centerPos,5,1);
                title('eye position across run');
            end
            
            % chop up by condition
            
            conds = [trial.cond]; conds(conds>0) = [condition(conds(conds>0)).stim]; conds = conds+1;
            c = struct('name',{'blank' 'inverted' 'upright'},'samples',[]);
            onsets = [trial.onset trial(end).onset+4];
            for n = 1:length(conds)
                ind = findBetween(samples(:,1),onsets(n),onsets(n+1));
                c(conds(n)).samples = [c(conds(n)).samples;samples(ind,:)];
            end
            
            for n = 1:length(c)
                if runFigs
                    subplot(1,4,n+1);
                    eyeInSpace_scan(c(n).samples,ppd,fixRad,centerPos,5,1); % outputs a centered samples in dva
                    title(c(n).name);
                end
                
                % aggregate
                subjc(n).name = c(n).name;
                subjc(n).samples = [subjc(n).samples;c(n).samples];
                subjc(n).runInd = [subjc(n).runInd repmat(r,1,length(c(n).samples))];
                [~,m]=vpnlSessions('fixPRF',subjs{s});
                subjc(n).totalRuns = m;
            end
            
            if runFigs
                superTitle(etName);
                if saveFigs
                    niceSave([figDir 'run_figs/' subjs{s} '/'],['fixPRF_run' num2str(r) preText]);
                    close all;
                end
            end
            
        else
            fprintf('*** Missing %s...\n',etMat);
        end
    end
    if ~isempty(subjc(1).samples)
    % main figure - collapsed across runs
    niceFig([.1 .1 .8 .5]);
    for n = 1:length(c)
        subplot(1,4,n);
        subjc(n).convertedSamples = eyeInSpace_scan(subjc(n).samples,ppd,fixRad,centerPos,5,1); % outputs a centered samples in dva
        title(c(n).name);
        % calc cond-wise std
        subjc(n).stdX = nanstd([subjc(n).convertedSamples(:,2)]);
        subjc(n).stdY = nanstd([subjc(n).convertedSamples(:,3)]);
        subjc(n).stdXY = [subjc(n).stdX subjc(n).stdY];
    end
    subplot(1,4,4);
    niceGroupedBars(vertcat(subjc.stdXY),[],{subjc.name},{'X std' 'Y std'});
    axis square;
    superTitle(['fixPRF subj ' subjs{s} ', etdata from ' num2str(length(unique(subjc(n).runInd))) ' of ' num2str(subjc(n).totalRuns) ' runs']);
    
    if saveFigs
        niceSave([figDir 'concat_runs/'],['fixPRF_' subjs{s}]);
    end
    % save data from this subject
    fprintf('Saving %s...\n',[matDir subjs{s} '/fixPRF_cdata_' subjs{s} preText '.mat']);
    save([matDir subjs{s} '/fixPRF_cdata_' subjs{s} preText '.mat'], 'subjc','ppd','centerPos');
    else
        fprintf('*** Skipping %s...\n',subjs{s});
    end
end
%if onLaptop playSound; end
%end