% analysis script for initial import and reading of eyelink files
% updated from prfRec2 to fixPRF5 scan

%%% ET NAME: always in the format subj/fixPRF_runNum - will be parsed later
clear all; close all;

subjs = prfSubjs;
saveFig = 1;
restart = 0;

for s = 1:length(subjs)
    
    % set dirs & files
    exptDir = pwd;% [raid  'invPRF/eyetracking/'];
    
    for n = 1:10
        % input
        etName = sprintf('%s/fixPRF_%i',subjs{s},n);
        etFile = [exptDir 'edfs/' etName '.edf'];
        
        %output
        ascFile = [exptDir 'ascs/' etName '.asc'];
        %dataFile = [exptDir 'data/' dataName  '.mat'];
        matFile = [raid 'invPRF/eyetracking/mats/' etName '_preprocessed'];
        figDir =  [exptDir 'figures/'];
        
        if ~exist([matFile '.mat']) || restart == 1
            checkDir(dirOf(matFile));
            if exist(etFile)
                fprintf('Starting %s...\n',etFile);
                % read asc file if necessary - this should be done at the end of the
                % eyetracking session
                if ~exist(ascFile,'file') edfRead([exptDir 'edfs/' etName '.edf'],[exptDir 'ascs/' subjs{s} '/']); end
                
                % parse ASC file to event & sample info
                %if ~exist(matFile,'file') || doParse [trial,info] = ascParse(ascFile,dataFile,matFile); else load(matFile); end
                %info
                
                % plot trials & condition averages
                %eyetrackPlots(exptDir,etName,dataName,0,[1 0 1]);
                
                % make some figures without parsing into trials yet
                sampleFile = [dirOf(ascFile) 'fixPRF_' num2str(n) '_samples.asc'];
                %niceFig(plt.figSize,18); niceGCA;
                [samples, startTime] = ascSampleRead(sampleFile); % samples = time, x, y pos in pixels
                % rate in core asc file == 250 samples/sec
                rate = 250;
                exptSec = 272;
                % GAZE_COORDS 0.00 0.00 1920.00 1080.00
                
                % if loadScanInfo % for now, just do this once per subject, since parsing these guys by trials is going to be...ugh
                % scanFile = [raid 'invPRF/fixPRF/' vpnlSessions('fixPRF',subjs{s}) '/Stimuli/output/fixPRF_' num2str(n) '.mat'];
                %     load(scanFile);
                %     ppd = scan.ppd;
                % loadScanInfo = 0;
                % end
                
                % eyeInTime.m function does this much more sanely
                
                if length(samples)/rate < exptSec
                    fprintf('Not enough samples, skipping...\n')
                else
                
                % set time to zero
                nativeTime = samples(:,1); % for indexing later
                samples(:,1) = (samples(:,1)-samples(1,1))/1000;
                dva = 0;
                
                plotText = {'X position ' 'Y position '}; if dva strcat(plotText,'(dva)'); else strcat(plotText,'(pix)'); end
                
                h1= niceFig([.1 .2 .9 .8]);
                
                for p = 1:2
                    % make plot
                    set(gca,'YDir','reverse'); subplot(4,2,p);
                    plot(samples(:,1),samples(:,p+1),'LineWidth',2,'Color',condColors(p));
                    xlim([0 samples(end,1)]);
                    title({plotText{p}; num2str(length(samples)/rate)}); xlabel('Trial Time (s)'); ylabel(plotText{p});
                end
                
                
                for p = 1:2
                    % make plot
                    set(gca,'YDir','reverse'); subplot(4,2,p+2);
                    processed = removeBlinks(samples);
                    plot(processed(:,1),processed(:,p+1),'LineWidth',2,'Color',condColors(p));
                    xlim([0 processed(end,1)]);
                    title('blinks removed'); xlabel('Trial Time (s)'); ylabel(plotText{p});
                end
                
                for p = 1:2
                    % make plot
                    set(gca,'YDir','reverse'); subplot(4,2,p+4);
                    processed = removeSpikes(processed);
                    plot(processed(:,1),processed(:,p+1),'LineWidth',2,'Color',condColors(p));
                    xlim([0 processed(end,1)]);
                    title('spikes removed'); xlabel('Trial Time (s)'); ylabel(plotText{p});
                end
                
                %%% trim to 272 - first randomly, then take input from user to
                %%% adjust
                startTrim = round(rand*(length(samples)/rate-exptSec)*rate);
                
                %aa = processed(startTrim:endTrim,:);
                %length(aa)/rate;
                
                h2= niceFig([.4 .3 .9 .7]);
                while 1
                    
                    endTrim = (startTrim+exptSec*rate)-1; % both of these are in ind of sample, not time units
                    
                    if endTrim > length(processed)
                        endTrim = length(processed);
                        startTrim = (endTrim-exptSec*rate)+1;
                        fprintf('>>> Selected beyond end of samples - autocorrecting.\n');
                    end
                    
                    for p = 1:2
                        % make plot
                        set(gca,'YDir','reverse'); subplot(2,1,p);
                        
                        plot(processed(:,1),processed(:,p+1),'LineWidth',2,'Color',condColors(p));
                        vline(startTrim/rate); vline(endTrim/rate);
                        xlim([0 processed(end,1)]);
                        title('proposed trim'); xlabel('Trial Time (s)'); ylabel(plotText{p});
                    end
                    figure(h2);
                    
                    trimOk = input('Trim ok? 1 = yes, 2 = rand, 3 = manual: ');
                    if trimOk==1 break;
                    elseif trimOk == 2
                        startTrim = round(rand*(length(samples)/rate-exptSec)*rate);
                    elseif trimOk == 3
                        [startTrim,~] = ginput(1);
                        startTrim = startTrim*rate;
                        
                    end
                end
                
                figure(h1);
                for p = 1:2
                    % make plot
                    set(gca,'YDir','reverse'); subplot(4,2,p+6);
                    
                    plot(processed(:,1),processed(:,p+1),'LineWidth',2,'Color',condColors(p));
                    vline(startTrim/rate); vline(endTrim/rate);
                    xlim([0 processed(end,1)]);
                    title('trim'); xlabel('Trial Time (s)'); ylabel(plotText{p});
                end
                
                fprintf('startTrim = %.2f, endTrim = %.2f.\n',startTrim,endTrim);
                samples = processed(startTrim:endTrim,:);
                
                if ~isequal(length(samples)/rate,exptSec)
                    error('Something is up with this trim length!'); end
                
                saveThis = input(['Save ' matFile '? 1/2: ']);
                if saveThis == 1
                    save(matFile,'samples','rate','sampleFile','startTrim','endTrim');
                end
                
                superTitle([etName],12,.05);
                if saveFig
                    niceSave([figDir 'preprocessing/' subjs{s} '/'],['fixPRF_' num2str(n) '_preprocessing']);
                    close all;
                end
                
                
                end
            else
                fprintf('*** Missing %s...\n',etFile);
            end
        else
            fprintf('Skipping %s, already done...\n',etFile);
        end
        
        
        %end
        %cd(exptDir);
    end
    %if onLaptop playSound; end
end