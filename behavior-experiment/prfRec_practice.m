function prfRec_practice(s,taskBack,face)
%clear all; close all;
% version 2: provides more detailed feedback to subject
%
%%% inputs:
% subject: string ID of subject (default to 'test')
% runNum: double, scan number (default to 999)
%
% adapted from scanner 2-back experiment to flexibly test behavior for
% upright/inverted faces at different retinal locations
% options: n-back task
% options: include all viewpoints or just front-view faces
% position is randomized between blocks
% inv/upright is randomized within blocks

%taskBack = 2;
expt.taskBack = taskBack;  % double, n-back recognition task

if expt.taskBack == 1
    expt.front = 0;     % binary, include only front-view faces or all views
else expt.front = 1;  end

expt.eyetrack = 1;  % binary, whether we are eyetracking this session
expt.plotGaze = 1;  % binary, whether we plot the current eye position on the screen
expt.driftCorr = 1; % binary, whether we use manual drift correction in the plotting

input('Hit enter to proceed.');
Screen('Preference', 'SkipSyncTests', 1);

%%% grab keyboard indices
[external, laptop] = deviceNums;

% to ensure that we know/use the correct monitor resolution, run
% ResolutionTest to see available monitor settings, and set it
expt.res = [1024 768 85];%[800,600,85];

%%%% set-up rand
rand('twister', sum(100*clock));
expt.rand = rand;

%%%% files and things
expt.root = pwd;
expt.date = datestr(now,30);
expt.edfName = ['prac'];

%%%% images used
% faceViews.mat contains struct of face images in format face(n).view(m) m=7 is the front-view;
if ~exist('face','var')
    fprintf('Loading stimuli...'); load('stims/faceViews.mat'); fprintf('Done.\n'); end
ims.numIDs = length(face);

%%%% timing
params.trialLength = 6;                 % in seconds, length of trials. number of faces is determined from this, timeBase & imSet
params.timeBase = .2;                   % in seconds
params.imSet = [1 1 1 1 0];             % timing (on/off) of each face in units of timeBase;
params.ITI = 1;                         % break between trials, in seconds
params.blockTrials = 6;                 % number of trials in each block (/2 orientations). must be even
params.locReps = 2;                     % number of blocks at each loc in this session

%%%% layout
params.faceSizeDeg = 3.2;               % size of each face (diam)
params.fixRadDeg =  .15;                % in degrees, the size of the biggest white dot in the fixation
params.outerFixPixels = 1;              % in pixels, the thickness of the cue ring. now draws *inside* fixRadDeg, so that the cue is not made bigger by increasing this param
params.fixAlpha = .25;                  % transparency of fixation point, exclusive of outer ring and middle point
params.midColor = [0 0 0];              % color of middle dot

%%%% experimental manipulation - cue locs
params.locsDeg = [[0 0];[-1 1];[1 -1]];         % in degrees, the locations that we will cue
params.stimNames = {'Inv-'; 'Upr-'};            % interpretable names of our conditions
params.locNames = {'Center'; 'Lower'; 'Upper'}; % interpretable names of our conditions

%%%% screen
params.textColor = [255 255 255];

%%% SET THESE IN TESTING ROOM
expt.screenWidth = 38.5;                % in cm; % et room = 38.5, laptop=33, CNI = 104cm at both 16ch and 32ch
expt.viewingDist = 54;                  % in cm; % et room = 54; laptop= ~40 (SP)

%%%% task
params.task.back = expt.taskBack;        % 1 or 2-back
params.task.text = {'Task: 1-Back on Faces.\n Press the Space Bar to indicate immediate repeats of the same IDENTITY, across different viewpoints.\n',...
    'Task: 2-Back on Faces.\n Press the Space Bar to indicate 2-back repeats of the same IDENTITY; all faces will be front-view.\n'};          % text for subject
params.task.prob = .2;                  % proportion of faces within trial that are targets
params.task.window = 5;                 % number of flips (timeBase units) in which we'll count a correct response

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               trial randomization                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% for easier indexing
numStim = 2;            % 2 stim types - inverted, normal
numLocs = length(params.locsDeg);
numConds = numLocs * numStim;

% struct of experimental conditions
condition = struct('pos',num2cell(repmat([1:numLocs],1,numStim)),...
    'stim',num2cell(Expand([1:numStim],numConds/numStim,1)),... % inv, normal
    'name',cellstr(strcat(Expand(params.stimNames,1,numConds/numStim),...
    strjust(strcat(repmat(params.locNames,numStim,1)),'left'))'));

% struct of trial specifics
% location in shuffled across blocks, stim is shuffled within blocks
block = struct('loc',num2cell(Shuffle(repmat([1:numLocs],1,params.locReps))));

for n = 1:length(block) % determine trial order etc
    
    % randomize stim/trial order
    trialOrder = Shuffle(repmat([1:numStim],1,params.blockTrials/numStim));
    
    % number of faces per trial
    trialSet = params.trialLength/(length(params.imSet)*params.timeBase);
    if ~mod(trialSet,1) == 0 % check if it's an integer
        error('trialLength/timeBase numbers don''t match up!');
    else trialSet = [1:trialSet];
    end
    
    % trial specifics
    block(n).trial = struct('stim',num2cell(trialOrder),'pos',num2cell(repmat(block(n).loc,1,params.blockTrials)),...
        'numTargs',[0],'targInd',zeros(1,length(trialSet)),'pressFlip',[]);
    
    for t = 1:params.blockTrials
        
        % assign condition number, for easier indexing later
        for c = 1:length(condition)
            if  block(n).trial(t).pos == condition(c).pos &&  block(n).trial(t).stim == condition(c).stim
                block(n).trial(t).cond = c; end
        end
        
        targPoss = [];
        % determine number of targets in this trial probibalistically
        block(n).trial(t).numTargs = length(find(rand(1,length(trialSet))<params.task.prob));
        if block(n).trial(t).numTargs == 0 block(n).trial(t).numTargs = 1; end
        
        % set targets
        while sum(block(n).trial(t).targInd) < block(n).trial(t).numTargs
            if isempty(targPoss) % reset this process if we run into an impossible targPoss situation
                
                if sum(block(n).trial(t).targInd)>0 % reset number of targets for good measure - prevents some weird loops
                    block(n).trial(t).numTargs = length(find(rand(1,length(trialSet))<params.task.prob));
                    if block(n).trial(t).numTargs == 0 block(n).trial(t).numTargs = 1; end
                end
                
                targPoss = find(trialSet); targPoss(1:params.task.back)=[]; targPoss(end) = [];
                block(n).trial(t).targInd = zeros(1,length(trialSet));
            end
            
            targInd = datasample(targPoss,1);
            block(n).trial(t).targInd(targInd) = 1; % indexing when targets occur
            
            % edit targPoss so we don't end up with chains of targets
            [~, ii, ~] = intersect(targPoss,targInd-params.task.back:targInd+params.task.back);%targPoss==targInd-params.task.back)=[]; targPoss(targPoss==targInd+params.task.back)=[];
            targPoss(ii) = [];
        end
        
        % set face IDs and views
        block(n).trial(t).IDs = datasample(1:ims.numIDs,length(trialSet),'Replace',false);
        if expt.front block(n).trial(t).views = repmat(7,1,length(trialSet));
        else block(n).trial(t).views = datasample(1:7,length(trialSet),'Replace',true); end
        
        for targInd = find(block(n).trial(t).targInd)
            block(n).trial(t).IDs(targInd) = block(n).trial(t).IDs(targInd-params.task.back); % set target equal to 2-back face
            % if this is a 1-back task & all views are used, don't repeat views consecutively
            if expt.taskBack && ~expt.front && block(n).trial(t).views(targInd) == block(n).trial(t).views(targInd-params.task.back)
                allViews = [1:7]; allViews(block(n).trial(t).views(targInd))=[];
                block(n).trial(t).views(targInd-params.task.back) = Sample(allViews);
            end
        end
        
        % indexing stimulus properties for each screenflip within the trial
        block(n).trial(t).flip.IDs = Expand(block(n).trial(t).IDs,length(params.imSet),1).*...
            repmat(params.imSet,1,length(trialSet)); % face ID - zeros are blanks
        block(n).trial(t).flip.views = Expand(block(n).trial(t).views,length(params.imSet),1).*...
            repmat(params.imSet,1,length(trialSet)); % fave view
        block(n).trial(t).flip.targ = Expand(block(n).trial(t).targInd,length(params.imSet),1); % blank periods are still valid target IDs for this trial
        block(n).trial(t).flip.im = Expand(trialSet,length(params.imSet),1); % which # face within the trial (for indexing)
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
expt.time = (params.blockTrials*params.locReps*numLocs)*(params.trialLength+params.ITI);   % approx time in seconds (not including inter-block breaks)
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
    [~,el,eyeInit] = initEyelink(expt.edfName);
    
    %%% initial window - wait for experimenter to advance
    [win,rect] = screenInit(params.backgroundColor,18,1); % launch screen at max(screens), enabling alpha channel
    DrawFormattedText(win,['Please keep your head still and in the chinrest during these practice trials.\n This ensures stable eyetracking.\n'...
        'Press Space to proceed.'], 'center','center',params.textColor);
    if expt.plotGaze gazeX = []; gazeY = []; end
end
Screen(win, 'Flip', 0);
KbQueueTriggerWait(external,'Space');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       experiment                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% instruction screen: task directions
DrawFormattedText(win,[params.task.text{params.task.back} 'Press Space to proceed.'], 'center','center',params.textColor);
Screen(win, 'Flip', 0);
KbQueueTriggerWait(external,'Space');


ptbDotWait(win,3)

expt.startTime = GetSecs;
tr = 1;
%%%%%%% START task TASK/FLIPPING
for n = 1:length(block)
    
    for t = 1:params.blockTrials
        gaze(tr).x = []; gaze(tr).y = [];
        %***!!!***!!!***!!!***!!!***!!!***!!!%
        % start trial on eyetracker
        if expt.eyetrack
            gaze(tr).x = []; gaze(tr).y = [];
            if  expt.driftCorr
                duration = 3; % in seconds, how long we record the eye (minus a 1s buffer that lets the subject fixate)
                [gaze(tr).corrX,gaze(tr).corrY, gaze(tr).SDs] = manualDriftCorrect(0,el,win,xc,yc,duration);
            else gaze(tr).corrX = 0; gaze(tr).corrY = 0;
            end
            eyelinkStartTrial('PRACTRIALSTART'); end
        %***!!!***!!!***!!!***!!!***!!!***!!!%
        
        for f = 1:(length(trialFlip)-1)
            
            if block(n).trial(t).flip.IDs(f) > 0 % not a blank period - draw a face
                fIm = face(block(n).trial(t).flip.IDs(f)).view{block(n).trial(t).flip.views(f)};
                if block(n).trial(t).stim == 1 % invert face
                    fIm = flipud(fIm); end
                fTex = Screen('MakeTexture',win,fIm);
                
                Screen('DrawTexture', win, fTex,[],expt.rects(block(n).trial(t).pos,:));
                Screen('Close', fTex);
            end
            
            % fixation
            bullseyeFix(win,xc, yc, expt.fixRad, params.outerFixPixels, params.midColor, params.fixAlpha,1);
            
            %***!!!***!!!***!!!***!!!***!!!***!!!%
            % show subject's current eye position
            if expt.plotGaze [gaze(tr).x ,gaze(tr).y ] = eyelinkGetPosition(el,gaze(tr).x ,gaze(tr).y);
            dots = 20; % length of samples that we will use for drawing
            if length(gaze(tr).x) > dots
                Screen('DrawDots',win,[gaze(tr).x(end-dots:end)+gaze(tr).corrX;gaze(tr).y(end-dots:end)+gaze(tr).corrY],2,[200 0 0]);
            end
            end
            %***!!!***!!!***!!!***!!!***!!!***!!!%
            
            %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if f == 1 % first flip of this trial
                if t == 1 % first trial of this block
                    [VBLT block(n).trial(t).startTime FlipT missed] = Screen(win, 'Flip', 0);
                else [VBLT block(n).trial(t).startTime FlipT missed] = Screen(win, 'Flip',ITIstart+params.ITI-slack);
                end
                block(n).trial(t).flipTime(f) = block(n).trial(t).startTime; % mark trial start
            else [VBLT block(n).trial(t).flipTime(f) FlipT missed] = Screen(win, 'Flip', block(n).trial(t).startTime + trialFlip(f) - slack);end % mark all flips
            
             %***!!!***!!!***!!!***!!!***!!!***!!!%
             % continue to monitor subject's eye position
            if expt.plotGaze
                while GetSecs < block(n).trial(t).startTime + params.timeBase - 2*slack % within this flip
                    [gaze(tr).x,gaze(tr).y] = eyelinkGetPosition(el,gaze(tr).x,gaze(tr).y); end
            end
            %***!!!***!!!***!!!***!!!***!!!***!!!%
             
            % listen for response
            [pressed, firstPress]= KbQueueCheck(external);
            if pressed == 1 block(n).trial(t).pressFlip = [block(n).trial(t).pressFlip f];
                KbQueueFlush();
            end
            
        end % end of trial
        
        % fixation
        bullseyeFix(win,xc, yc, expt.fixRad, params.outerFixPixels, params.midColor, params.fixAlpha,1);
        
        %%%% to show the very last flip screen for its 100ms
        [VBLT block(n).trial(t).flipTime(end+1) FlipT missed] = Screen(win, 'Flip',block(n).trial(t).startTime + trialFlip(end)  - slack);
        
        %***!!!***!!!***!!!***!!!***!!!***!!!%
        if expt.eyetrack eyelinkEndTrial('PRACTRIAL_END'); end
        %***!!!***!!!***!!!***!!!***!!!***!!!%
        
        % fixation & ITI
        if block(n).trial(t).flip.IDs(end) > 0 % ITI shows fixation, but it should already be drawn
            bullseyeFix(win,xc, yc, expt.fixRad, params.outerFixPixels, params.midColor, params.fixAlpha,1);
            
            [VBLT ITIstart FlipT missed] = Screen(win, 'Flip', 0);
        else ITIstart = GetSecs;
        end
        tr = tr+1;
    end
end

%%%%%%%%%%%%%%%%%%
% done! wrap up  %
%%%%%%%%%%%%%%%%%%

[VBLT ITIstart FlipT missed] = Screen(win, 'Flip', 0);
expt.runTime = GetSecs - expt.startTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% performance
perf = struct('condName',{condition.name},'numTargs',0,'hits',[],'FAs',[],'hitRate',[],'FArate',[]);

% aggregate hits and false alarms by condition/block
% good god, this is hack-y.

for n = 1:length(block)
    for t=1:params.blockTrials
        hits = []; FAs = [];
        for p = block(n).trial(t).pressFlip
            targetRange = p-params.task.window:p-1;
            try lastTarg = find(block(n).trial(t).flip.targ(targetRange),1,'last');
            catch
                lastTarg = [];
            end
            if ~isempty(lastTarg)
                hits = [hits block(n).trial(t).flip.im(targetRange(lastTarg))]; % this is a hit
            else
                FAs = [FAs p]; % this is a false alarm
            end
        end
        perf(block(n).trial(t).cond).hits = [perf(block(n).trial(t).cond).hits unique(hits)]; % only one press per target counts
        perf(block(n).trial(t).cond).FAs = [perf(block(n).trial(t).cond).FAs FAs];
        perf(block(n).trial(t).cond).numTargs = perf(block(n).trial(t).cond).numTargs + block(n).trial(t).numTargs;
    end
end

for c = 1:length(perf)
    perf(c).hitRate = length([perf(c).hits])/perf(c).numTargs;
    perf(c).FArate = length([perf(c).FAs])/(params.locReps*(params.blockTrials/numStim)*length(trialSet)-perf(c).numTargs);
    
    %%% do we need to apply correction to this dprime rate?
    % correct for 0 false alarms by setting it to 0.5 false alarms, and 1
    % hit rate by setting it to 0.5 errors
    txt = 'uncorrected';
    if perf(c).FArate == 0 FAfix = .5/(params.locReps*(params.blockTrials/numStim)*length(trialSet)-perf(c).numTargs); txt = '*corrected'; else FAfix = perf(c).FArate; end
    if perf(c).hitRate == 1 HRfix = 1-(.5/perf(c).numTargs); txt = '*corrected'; else HRfix = perf(c).hitRate; end
    
    perf(c).dprime = norminv(HRfix)-norminv(FAfix);
    
    fprintf('%s Condition: %0.2f hit rate, %0.2f false alarm rate, %s d'' = %0.2f.\n',perf(c).condName,perf(c).hitRate,perf(c).FArate,txt,perf(c).dprime);
end


if expt.eyetrack
    eyelinkClose;
    allGaze = gaze;
    %%% give feedback to participant
    
    for n = 1:length(gaze)
        gaze(n).xSD = nanstd([gaze(n).x])/expt.ppd;
        gaze(n).ySD = nanstd([gaze(n).y])/expt.ppd;
        
        if expt.driftCorr
            gaze(n).x = gaze(n).x+gaze(n).corrX;
            gaze(n).y = gaze(n).y+gaze(n).corrY;
        end
    end
    
    %%% eyetracking feedback
    feedbackTr = round(length(gaze)/2):length(gaze); % base feedback on second half of practice trials
    
    for t = feedbackTr
    Screen('DrawDots',win,[[gaze(t).x];[gaze(t).y]],2,condColors(t,1)*255);
    end
    bullseyeFix(win,xc, yc, expt.fixRad, params.outerFixPixels, params.midColor, params.fixAlpha,1);% Screen(win, 'Flip',0);
    
    medXdev = nanmedian([gaze(feedbackTr).xSD]); medYdev = nanmedian([gaze(feedbackTr).ySD]);
    if medXdev < .5 perfText = 'Good!';
    elseif medXdev < 1 perfText = 'Average.';
    else perfText = 'Fair. Recommend more practice trials.'; end
    
    perfText = sprintf('Your eye position over the second half of the practice session is shown below.\nAcross these trials, your median deviation was %.2fdva in the X direction\nand %.2fdva in the Y direction.\nThis fixation performance is %s\n',...
        medXdev,medYdev,perfText);
    
    fprintf('medXdev = %.2f, medYdev = %.2f. %s\n', medXdev, medYdev,perfText);
    
    % drift-corr fixation performance - shown only to experimenter
    fprintf('median SDs during driftCorr = %.2fdva\n', nanmedian([gaze(feedbackTr).SDs]./expt.ppd))
    
    DrawFormattedText(win,[perfText '\nPress Space to proceed.'], 'center',yc *.25 ,params.textColor);
    Screen(win, 'Flip', 0);fixIm=Screen('GetImage', win); imwrite(fixIm,['pracFix/' s '_practiceFix.png'],'png');
    
    KbQueueTriggerWait(external,'Space');
    
    %%% behavioral feedback
    feedbackCond = 4; % make sure this is upright-center
    
    if perf(feedbackCond).hitRate > .7 perfText = 'Great!';
    elseif perf(feedbackCond).hitRate > .4 perfText = 'Good.';
    else perfText = 'Fair. Recommend more practice trials.'; end
    fprintf('%s hit rate: = %.2f%%, num FAs = %d. %s\n',perf(feedbackCond).condName,...
        perf(feedbackCond).hitRate*100, length(perf(feedbackCond).FAs),perfText);
    
    perfText = ['You detected ' (num2str(perf(feedbackCond).hitRate*100)) '% of the repeats and made ' num2str(length(perf(feedbackCond).FAs)) ' false-alarm keypresses\n'...
        'when the faces were upright and centrally presented.\nThis performance is ' perfText];
    DrawFormattedText(win,[perfText '\nPress Space to proceed.'], 'center','center',params.textColor);
    Screen(win, 'Flip', 0);
    KbQueueTriggerWait(external,'Space');
    
end

KbQueueRelease();
fclose all;
sca;
%end
playSound;



