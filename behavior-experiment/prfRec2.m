function [params,expt,block,perf,eyeInit] = prfRec2(initials,face)
% version 2: 2AFC trials with paired faces presentation 
%clear all; close all;
sca;

% check for path to sputils
if ~exist('alphabet.m','file'); addpath(genpath('~/Experiments/Sonia/matlab/utils'));end

%%% inputs:
% subject: initials (default to '00')
% face: mat of face images - to save us some loading time, optionally
%
% adapted from scanner 2-back experiment to flexibly test behavior for
% upright/inverted faces at different retinal locations
% options: which view points are we using?
% options: number and timing of faces: currently using 3 400ms faces +
% 400ms cue, 800ms delay

expt.views = [7 1 4];  % which viewpoints from faceViews.mat will be used?
% 7 = front, order = [6 5 4 7 1 2 3]; see faceViews.fig

%%% settings for this session's eyetracking
expt.eyetrack = 0;  % binary, whether we are eyetracking this session
expt.driftCorr = 1; % binary, whether we do drift correction every block
expt.plotGaze = 0;  % binary, whether we plot the current eye position on the screen
expt.feedback = 0;  % giving feedback on the same/different task

input('Hit enter to proceed.');
Screen('Preference', 'SkipSyncTests', 1);

%%% grab keyboard indices
[external, laptop] = deviceNums;

% to ensure that we know/use the correct monitor resolution, run
% ResolutionTest to see available monitor settings, and set it
if onLaptop expt.res = [];
else expt.res = [1024 768 85];end   %[800,600,85];

%%% default vars
if ~exist('initials','var')
    expt.subj = '00';
else expt.subj = initials;
end

%%%% set-up rand
rand('twister', sum(100*clock));
expt.rand = rand;

%%%% files and things
expt.root = pwd;
expt.date = datestr(now,30);
expt.edfName = ['prfRec2_' expt.subj '_' datestr(now,'yymmdd')];
expt.soundOn = 1;

%%%% images used
% faceViews.mat contains struct of face images in format face(n).view(m) m=7 is the front-view;
if ~exist('face','var')
    fprintf('Loading stimuli...'); load('stims/faceViews.mat'); fprintf('Done.\n'); end
load('stims/facePairs.mat'); % struct pairs{n} = two matched faces
%load('stims/masks.mat'); % struct mask{n} = phase-scrambled masking images
load('stims/feedbackSounds.mat') % hit and miss sounds

%%%% timing
params.timeBase = .200;                 % in seconds (small so that we can draw eye position effectively)
params.faceDur = .400;                  % in seconds
params.numFaces = 3;                    % number of faces in encoding stream
params.cueDur = .4;                     % duration of cue face, in seconds
params.delay = .8;                      % in seconds
params.ITI = 2;                         % length of drift correction duration, too!! (in seconds, minus buffer specified in manualDriftCorrection.m
params.blockTrials = 20;                % number of trials in each block (/2 orientations). must be /4 for equal numbers of hits and misses
params.locReps = 6;                     % number of blocks at each loc in this session
params.condTrials = params.locReps*params.blockTrials/2;    % number of trials in each condition
expt.task = sprintf('Stigliani-based WM task with %d faces in stream presented at %.0fms, %d trials per condition)',...
    params.numFaces,params.faceDur*1000,params.condTrials);
trialSet = [Expand([1:params.numFaces],params.faceDur/params.timeBase,1)...
    zeros(1,params.delay/params.timeBase) repmat(params.numFaces+1,1,params.cueDur/params.timeBase)];

%%%% layout
params.faceSizeDeg = 3.2;               % size of each face (diam)
params.fixRadDeg =  .15;                % in degrees, the size of the biggest white dot in the fixation
params.outerFixPixels = 1;              % in pixels, the thickness of the cue ring. now draws *inside* fixRadDeg, so that the cue is not made bigger by increasing this param
params.fixAlpha = .25;                  % transparency of fixation point, exclusive of outer ring and middle point
params.midColor = [0 0 0];              % color of middle dot

%%%% experimental manipulation - cue locs
params.locsDeg = [[-.67 .79];[0 0];[+.67 -.79]];         % in degrees, the locations that we will cue
params.stimNames = {'Inv-'; 'Upr-'};                     % interpretable names of our conditions
params.locNames = {'LowerLeft'; 'Center'; 'UpperRight'}; % interpretable names of our conditions

%%%% screen
params.textColor = [255 255 255];

%%% SET THESE IN TESTING ROOM
expt.screenWidth = 38.5;                % in cm; % et room = 38.5, laptop=33, CNI = 104cm at both 16ch and 32ch
expt.viewingDist = 54;                  % in cm; % et room = 54; laptop= ~40 (SP)

%%%% task
params.task.keys = {'J' 'F'};
params.task.text = ['After each trial, report whether the stream of faces contained the target face.\nPress ' params.task.keys{1} ...
    ' = YES or ' params.task.keys{2} ' = NO.\n'];          % text for subject
params.task.keys = lower(params.task.keys);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               trial randomization                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% for easier indexing
numStim = 2;            % 2 stim types - inverted, normal
numLocs = length(params.locsDeg);
numConds = numLocs * numStim;

% approximate
fprintf('Experiment will take approx. %.2f mins.\n',...
    ((params.faceDur*params.numFaces+params.cueDur+params.delay+params.ITI)...   % single trial
    *params.blockTrials*params.locReps*length(params.locNames))/60);            % all the trials

% struct of experimental conditions
condition = struct('pos',num2cell(repmat([1:numLocs],1,numStim)),...
    'stim',num2cell(Expand([1:numStim],numConds/numStim,1)),... % inv, normal
    'name',cellstr(strcat(Expand(params.stimNames,1,numConds/numStim),...
    strjust(strcat(repmat(params.locNames,numStim,1)),'left'))'));

% counterbalancing trials in each block: half upright/inverted, half
% same/different
stim = repmat([1:numStim],1,params.blockTrials/numStim);
targ = Expand([0 1],length(stim)/2,1);

% location in shuffled across blocks, stim is shuffled within blocks
block = struct('loc',num2cell(Shuffle(repmat([1:numLocs],1,params.locReps))));

for n = 1:length(block) % determine trial order etc
    
    % randomize stim/trial order
    tRand = Shuffle(1:length(stim));
    
    % trial specifics
    block(n).trial = struct('stim',num2cell(stim(tRand)),'loc',num2cell(repmat(block(n).loc,1,params.blockTrials)),...
        'pair',num2cell(datasample([1:length(pair)],length(tRand),'Replace',false))...
        ,'IDs',[],'views',[],'targ',num2cell(targ(tRand)),'response',[],'correct',0,'corrX',0,'corrY',0);
    % in this version, IDs 3 is masking face
    
    for t = 1:params.blockTrials
        
        % assign condition number, for easier indexing later
        for c = 1:length(condition)
            if  block(n).trial(t).loc == condition(c).pos &&  block(n).trial(t).stim == condition(c).stim
                block(n).trial(t).cond = c; end
        end
        
        % set views
        if length(expt.views) > 1
            while 1 % don't repeat same view on very last encoding face
                block(n).trial(t).views = datasample(expt.views,max(trialSet)); 
            if diff(block(n).trial(t).views(end-1:end))>0
                break
            end
            end
            
        % if we're just using 1 view, this doesn't matter
        else block(n).trial(t).views = repmat(expt.views,1,max(trialSet)); end
        
        % set faceIDs
        targs = Shuffle(pair(block(n).trial(t).pair,:));
        foils = pair(1:end ~= block(n).trial(t).pair,:);% all pairs but the target
        
        block(n).trial(t).IDs = Shuffle([datasample(foils(:),params.numFaces-1,'Replace',false);targs(1)]);
        
        if block(n).trial(t).targ block(n).trial(t).IDs(end+1) = targs(1); % it target, last face is targs(1)
        else block(n).trial(t).IDs(end+1) = targs(2); end % else it's the paired face
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      open window                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if ~expt.eyetrack HideCursor; end
Priority(9);

%%%% open screen
sampleFace = face(1).view{1};
params.backgroundColor = sampleFace(1,1); % match color to our face images
[win,rect] = screenInit(params.backgroundColor,18,1,expt.res); % launch screen at max(screens), enabling alpha channel

%%%% timing optimization
expt.flipInt = Screen('GetFlipInterval',win);
slack = expt.flipInt/2;
expt.time = (params.blockTrials*params.locReps*numLocs)*(+params.ITI);   % approx time in seconds (not including inter-block breaks)
params.trialLength = params.timeBase*length(trialSet);
trialFlip = (0:params.timeBase:params.trialLength);

%%%% scale the stims for the screen
expt.ppd = pi* rect(3) / (atan(expt.screenWidth/expt.viewingDist/2)) / 360;
expt.faceSize = round(params.faceSizeDeg*expt.ppd);                 % in degrees, the size of our faces
expt.fixRad = round(params.fixRadDeg*expt.ppd);

%%% move xc and yc off screen center if needed (more relevant for the 7T scanner...)
xc = rect(3)/2;
yc = rect(4)/2;

%%% locations to probe
expt.centers = (params.locsDeg*expt.ppd)+repmat([xc yc],size(params.locsDeg,1),1);
expt.rects = CenterRectOnPoint([0 0 expt.faceSize expt.faceSize],expt.centers(:,1), expt.centers(:,2));

%%% make textures
for n = 1:length(block)
    for m = 1:length(block(n).trial)
        for f = 1:length(block(n).trial(t).IDs)
            fIm = face(block(n).trial(m).IDs(f)).view{block(n).trial(m).views(f)};
                if block(n).trial(m).stim == 1 % invert face
                    fIm = flipud(fIm); end
                block(n).trial(m).fTex{f} = Screen('MakeTexture',win,fIm); 
        end
    end
end
clear face;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           eyetracking init & calibration                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% start recording keys
KbQueueCreate(external);
KbQueueStart();

%***!!!***!!!***!!!***!!!***!!!***!!!%
if expt.eyetrack
    % initialize eyetracking, run calibration
    [~,el,eyeInit] = initEyelink('tmp');
    
    %%% initial window - wait for experimenter to advance
    [win,rect] = screenInit(params.backgroundColor,18,1); % launch screen at max(screens), enabling alpha channel
    DrawFormattedText(win,['Eyetracking calibration complete.\n. Please keep your head still and in the chinrest throughout the rest of the study.\n'...
        'Fixate the center point whenever it is on the screen.\nPress Space to proceed.'], 'center','center',params.textColor);
else
    eyeInit = 'No Tracking'; % for final output;
    DrawFormattedText(win,['Press Space to proceed.'], 'center','center',params.textColor);
end
Screen(win, 'Flip', 0);

KbQueueTriggerWait(external,'Space');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       experiment                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% instruction screen: task directions
DrawFormattedText(win,[params.task.text 'Press Space to proceed.'], 'center','center',params.textColor);
Screen(win, 'Flip', 0);
KbQueueTriggerWait(external,'Space');
% fixation

bullseyeFix(win,xc, yc, expt.fixRad, params.outerFixPixels, params.midColor, params.fixAlpha,1);
Screen(win, 'Flip', 0);

WaitSecs(1);

bullseyeFix(win,xc, yc, expt.fixRad, params.outerFixPixels, [0 167 0], params.fixAlpha,1);
Screen(win, 'Flip', 0);
WaitSecs(2);

expt.startTime = GetSecs;

%%%%%%% START task TASK/FLIPPING
for n = 1:length(block)
    
    if n>1
        WaitSecs(1);
        DrawFormattedText(win,['You have completed ' num2str(n-1) ' of ' num2str(params.locReps * numLocs) ' blocks.\n' ...
            'Take a break if you need it. \n\nPress Space to Advance, and don''t forget to FIXATE throughout each trial.']...
            ,'center', 'center',[0 0 0]);
        Screen(win, 'Flip', 0);
        fprintf('Subject has completed block %d of %d...\n',n-1,params.locReps * numLocs);
        KbQueueTriggerWait(external,'Space');
    end
    for t = 1:params.blockTrials
       
        %***!!!***!!!***!!!***!!!***!!!***!!!%
        % start trial on eyetracker
        if expt.eyetrack
            if  expt.driftCorr && t==1 % do this before the first trial
                [block(n).trial(t).corrX,block(n).trial(t).corrY,block(n).trial(t).SDs] = manualDriftCorrect(0,el,win,xc,yc,params.ITI-.5);
            end
            
            block(n).trial(t).gazeX = []; block(n).trial(t).gazeY = [];
            eyeInit.trialMessage = 'BLOCK %d, TRIAL %d, CONDITION %d';
            eyelinkStartTrial(sprintf(eyeInit.trialMessage, n, t, block(n).trial(t).cond));
        end
        
        %***!!!***!!!***!!!***!!!***!!!***!!!%
        for f = 1:(length(trialFlip)-1)
            
            if trialSet(f) > 0 % not a mask  draw a face
                fTex = block(n).trial(t).fTex{trialSet(f)};
                Screen('DrawTexture', win, fTex,[],expt.rects(block(n).trial(t).loc,:));
                if trialSet(f) == max(trialSet) % mark the cue
                    Screen('FrameOval',win,[0 167 0],expt.rects(block(n).trial(t).loc,:),1,1); end
            else
            end
            
            % fixation
            bullseyeFix(win,xc, yc, expt.fixRad, params.outerFixPixels, params.midColor, params.fixAlpha,1);
            
            %***!!!***!!!***!!!***!!!***!!!***!!!%
            % show subject's current eye position
            if expt.plotGaze && expt.eyetrack
                [block(n).trial(t).gazeX,block(n).trial(t).gazeY] = eyelinkGetPosition(el,block(n).trial(t).gazeX,block(n).trial(t).gazeY,win,block(n).trial(t).corrX,block(n).trial(t).corrY);
            elseif expt.eyetrack % gather gaze data but don't plot it
                [block(n).trial(t).gazeX,block(n).trial(t).gazeY] = eyelinkGetPosition(el,block(n).trial(t).gazeX,block(n).trial(t).gazeY);
            end
            %***!!!***!!!***!!!***!!!***!!!***!!!%
            
            %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if f == 1 % first flip of this trial
                
                if t == 1 % first trial of this block - start now
                    [~,block(n).trial(t).startTime,~,~] = Screen(win, 'Flip', 0);
                    
                else
                    % starting trial relative to previous ITIstart
                    [~,block(n).trial(t).startTime,~,~] = Screen(win, 'Flip',ITIstart+params.ITI-slack);
                end
                block(n).trial(t).flipTime(f) = block(n).trial(t).startTime; % mark trial start
            else % not first flip of the trial
                [~,block(n).trial(t).flipTime(f),~,~] = Screen(win, 'Flip', block(n).trial(t).startTime + trialFlip(f) - slack);
            end
           
            %***!!!***!!!***!!!***!!!***!!!***!!!%
            % get subject's current eye position
            if expt.eyetrack
                [block(n).trial(t).gazeX,block(n).trial(t).gazeY] = eyelinkGetPosition(el,block(n).trial(t).gazeX,block(n).trial(t).gazeY);
            end
            %***!!!***!!!***!!!***!!!***!!!***!!!%
            
        end % end of trial
        
        
        % fixation
        bullseyeFix(win,xc, yc, expt.fixRad, params.outerFixPixels, params.midColor, params.fixAlpha,1);
        
        %%%% to show the very last flip screen for its 200ms
        [VBLT block(n).trial(t).flipTime(end+1) FlipT missed] = Screen(win, 'Flip',block(n).trial(t).startTime + trialFlip(end)  - slack);

        % fixation & ITI
        ITIstart = GetSecs; buffer = .5; % time in the ITI that we'll use for other things, like response checking

        %***!!!***!!!***!!!***!!!***!!!***!!!%
        if expt.eyetrack eyelinkEndTrial('TRIAL_END'); end
        if expt.driftCorr && t < length(block(n).trial) && expt.eyetrack
           [block(n).trial(t+1).corrX,block(n).trial(t+1).corrY,block(n).trial(t+1).SDs] = manualDriftCorrect(0,el,win,xc,yc,params.ITI-buffer);
        else WaitSecs(params.ITI-buffer); end
       
        %***!!!***!!!***!!!***!!!***!!!***!!!%
        % check response
                [pressed, firstPress]= KbQueueCheck(external);
                if pressed == 1
                    if ~isempty(cellNum(KbName(firstPress),params.task.keys)) % if pressed == one of the relevant keys
                        block(n).trial(t).response = abs(cellNum(KbName(firstPress),params.task.keys)-2); % recode to 1 = 1, 2 = 0
                        if block(n).trial(t).response == block(n).trial(t).targ block(n).trial(t).correct = 1; end
                        % auditory feedback
                        if expt.feedback
                            if block(n).trial(t).correct sound(hSound); else sound(mSound); end
                        end
                    end
                end
                KbQueueFlush();
    end    

        % interim save
        checkDir([pwd '/data']);
        if expt.eyetrack
            eval(['save data/prfRec2_' expt.subj '.mat params expt condition block eyeInit el']);
        else
            eval(['save data/prfRec2_' expt.subj '.mat params expt condition block']);
        end
end

%%%%%%%%%%%%%%%%%%
% done! wrap up  %
%%%%%%%%%%%%%%%%%%

[VBLT ITIstart FlipT missed] = Screen(win, 'Flip', 0);
[pressed, firstPress]= KbQueueCheck(external);
if pressed == 1
    KbName(firstPress)
    if ~isempty(cellNum(KbName(firstPress),params.task.keys)) % if pressed == one of the relevant keys
        block(n).trial(t).response = abs(cellNum(KbName(firstPress),params.task.keys)-2); % recode to 1 = 1, 2 = 0
        if block(n).trial(t).response == block(n).trial(t).targ block(n).trial(t).correct = 1; end
        if expt.feedback if block(n).trial(t).correct playSound; end
        end
    end
end
KbQueueFlush();
KbQueueRelease();

expt.runTime = GetSecs - expt.startTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DrawFormattedText(win,['You have completed the experiment!']...
    ,'center', 'center',[0 0 0]);
WaitSecs(1);
if expt.soundOn playSound; end

Screen(win, 'Flip', 0); t = GetSecs;

if expt.eyetrack
    eval(['save data/prfRec2_' expt.subj '_' expt.date '.mat params expt condition block perf eyeInit el']);
    eval(['save data/prfRec2_' expt.subj '.mat params expt condition block perf eyeInit el']);
else
    eval(['save data/prfRec2_' expt.subj '_' expt.date '.mat params expt condition block perf']);
    eval(['save data/prfRec2_' expt.subj '.mat params expt condition block perf']);
end

WaitSecs(5-(GetSecs-t)); % wait 5 seconds, exclusive of the time it took to save things

KbQueueRelease();

ShowCursor;
Screen('Close');
Screen('CloseAll');

% close graphics window, close data file and shut down tracker
if expt.eyetrack
    eyelinkGetEDF('tmp',expt.edfName,['~/Experiments/prfRec/edfs/']);
    eyelinkClose; end

sca;

%fclose all;
%clear all;
%end



