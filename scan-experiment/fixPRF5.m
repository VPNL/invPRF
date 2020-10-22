%function fixPRF5(subject,runNum,loadTrials,offset)
% assuming 1 rep of all conditions, expt: 191 (5TRs/10s of countdown)
% follows kay, weiner, grill-spector (2015) task experiment, presents faces
% either upright or inverted
% fixPRF: only fixation task
% version 5: trial-based fixation task on small set of letters

% clear all;
atScanner = 0;

% SP 8/2018: updates fixPRF to
% switch from 2-back to 1-back
% introduce flicker
% SP 9//13
% change flicker rate
% remove fixation grid
% make fixation non-transparent
% make fixation ring red
% also outputs performance plot for debugging
% ver 5: adds trial structure to RSVP string (rather than one
%
%%% inputs:
% subject: string ID of subject (default to 'test')
% runNum: double, scan number (default to 999)
% trials: optional, .mat file specifying the trial ordering to present -
%   useful if we need to match order etc across tasks/runs/subjects
% offset: optional, [x y] by which to offset screen center

%%% at the scanner:
% set trigger box to these options:
% change modes / manual configure / HHSC 1X5D / USB / HID NAR NO5

%%% for testing
%clear all;

scan.trigger = atScanner;   % set to zero to debug
scan.coil = 32 ; % 32 or 16ch coil at SNI
if scan.coil == 16 flipScreen = 1; else flipScreen = 0; end
input('Hit enter to proceed.');
Screen('Preference', 'SkipSyncTests', abs(atScanner-1));

[keyboardIndices, productNames, ~] = GetKeyboardIndices;
scanner = keyboardIndices(1);
if numel(keyboardIndices) == 1
    laptop = keyboardIndices(1); % don't differentiate scanner/laptop
else laptop = keyboardIndices(2); end


if ~exist('subject','var')
    scan.subj = 'test'; else scan.subj  = subject;
end

if ~exist('runNum','var')
    scan.runNum = 777; else scan.runNum  = runNum;
end

if ~exist('offset','var') || isempty(offset)
    scan.offset = [0 0]; else scan.offset = offset;
end

if ~exist('loadTrials','var')
    scan.loadTrials = []; else scan.loadTrials = loadTrials;
end

%%%% set-up rand
rand('twister', sum(100*clock));
scan.rand = rand;

%%%% files and things
scan.root = pwd; %'/Users/Sonia/Desktop/ObjectRF/';
scan.date = datestr(now,30);

%%%% images used
% faceIms.mat contains struct of face images in format face(n).view(m);
load('faceFront.mat');
ims.numIDs = 95;

%%%% timing
params.countDown = 10;                  % countdown at the beginning of the scan - 10
params.timeBase = .5;                   % in seconds - this is the letter length, faceID length
params.trialSet = [0 ones(1,7)];        % timing (on/off) of each trial in units of timeBase; 0 = cue onset
params.flickerBase = .1;                % timebase for intra-stimulus changes (on/off flicker)
params.flickerSet = [0 0 ones(1,3)];    % flicker is 2xflickerTime off 3xflickerTime on
params.initialNull = 4;                 % null trials at initial fixation - 4
params.finalNull = 4;                   % null trials at final fixation - 4
params.TRlength = 2;                    % in seconds
params.nullTrials = 10;                 % number of null trials randomly intermixed to estimate baseline activity

%%%% layout
params.grid = 5;                        % number of grid elements
params.gridSpaceDeg = 1.5;              % center to center spacing of the grid
params.faceSizeDeg = 3.2;               % size of each face (diam)
params.fixRadDeg =  .2;                 % in degrees, the size of the biggest white dot in the fixation
params.cuePix = 2;                      % in pixels, the thickness of the cue ring. now draws *inside* fixRadDeg, so that the cue is not made bigger by increasing this param

%%%% screen
%params.backgroundColor = [128 128 128]; % color - set later from sample image
params.textColor = [255 255 255];
scan.screenWidth = 104;                  % in cm; % laptop=27.5, office=43, CNI = 104cm at both 16ch and 32ch
if scan.coil == 32
    scan.viewingDist = 272;
else                                     % in cm; CNI = 270-273cm at 32ch, 265 at 16ch;
    scan.viewingDist = 265; end

%%%% tasks
params.taskBack = 1;                                                        % 1 or 2-back
params.task = [num2str(params.taskBack) '-Back on Letters'];                % for indexing
params.taskColor = [200 0 0];                                               % dark red
params.taskProb = .5;                                                       % proportion of trials with 2back ID targets
params.respWindow = 4*length(params.flickerSet);                            % number of flips in which we'll count a correct response
params.str = ['A' 'S' 'D' 'F' 'G' 'H' 'J' 'K'];% alphabet([],1);%           % entire alphabet, capitalized

%%%%%%%%%%%%%%%
% trial setup %
%%%%%%%%%%%%%%%

numStim = 2;            % 2 stim types - inverted, normal
numConds = params.grid * params.grid * numStim;

condition = struct('pos',num2cell(repmat([1:params.grid*params.grid],1,numStim)),...
    'stim',num2cell(Expand([1:numStim],numConds/numStim,1)),... % inv, normal
    'name',cellstr(strcat(Expand({'Inv-';'Norm-'},1,numConds/numStim),...
    strjust(num2str(repmat([1:params.grid*params.grid]',numStim,1)),'left'))'));

numTrials = numConds+params.initialNull+params.finalNull+params.nullTrials;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in this verion,targ 1 = RSVP targs, targ 2 = face targs
% RSVP targs happen on blank trials as well
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(scan.loadTrials)
    
    %randomly intersperses blank trials for baseline estim, pads with init/final fixation period
    trialOrder = [zeros(params.initialNull,1);...
        Shuffle([[1:length(condition)]';zeros(params.nullTrials,1)]); zeros(params.finalNull,1)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     % for debugging, a non-random trial order
    %     trialOrder = [zeros(params.initialNull,1);...
    %         [1:length(condition)]';zeros(params.nullTrials,1); zeros(params.finalNull,1)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    trial = struct('cond',num2cell(trialOrder),'onset',[],...
        'targ',[0 0],'targInd',[0 0]);
    
    % prep -  first and second face/letter can't be a target, by def of 2-back
    targPoss = find(params.trialSet); targPoss(1:params.taskBack)=[];
    
    % selects face identities, viewpoints, location of target if necessary
    for n = 1:numTrials
        trial(n).onset = (n-1)*(length(params.trialSet)*params.timeBase);
        
        %%%% face targets
        if trial(n).cond > 0 % if not blank
            if rand<params.taskProb trial(n).targ(2) = 1;end
            trial(n).IDs = datasample(1:ims.numIDs,length(params.trialSet),'Replace',false).*params.trialSet;
            if trial(n).targ(2) == 1 % determine index of face target
                trial(n).targInd(2)=datasample(targPoss,1); % which of our images in the trial will repeat the previous face
                trial(n).IDs(trial(n).targInd(2)) = trial(n).IDs(trial(n).targInd(2)-params.taskBack); % set target equal to 2-back face
            end
        else trial(n).IDs = zeros(1,length(params.trialSet));
        end
        %%%% letter targets - even if cond = 0
        if rand<params.taskProb trial(n).targ(1) = 1;end
        trial(n).rsvp = datasample(1:length(params.str),length(params.trialSet),'Replace',false).*params.trialSet; % 1 = blank
        if trial(n).targ(1) == 1 % determine index of RSVP target
            trial(n).targInd(1) = datasample(targPoss,1); % which of our images in the trial has the target dot
            trial(n).rsvp(trial(n).targInd(1)) = trial(n).rsvp(trial(n).targInd(1)-params.taskBack); % set target equal to 2-back letter
        end
        
    end
    
    save(['trialSets/' datestr(now,'mmdd') '_trials' num2str(scan.runNum) '.mat'],'trial');
else
    load(['trialSets/' scan.loadTrials]);
end

% indexing stimulus properties for each screenflip - now including flicker
% timing
flip.trials = [-1*Expand(fliplr([1:params.countDown]),1/params.flickerBase,1) ...
    Expand([1:numTrials],length(params.trialSet)*length(params.flickerSet),1)];
flip.targs = zeros(2,length(flip.trials));
flip.IDs = zeros(1,length(flip.trials)); flip.rsvp = flip.IDs; flip.flicker = flip.IDs;

% during expt, pos and task are grabbed from cond struct,
% ID/view/targetInd are grabbed from this flip struct

for n = 1:numTrials
    ind = find(flip.trials==n);
    flip.flicker(ind) = repmat(params.flickerSet,1,length(params.trialSet));
    flip.IDs(ind) = Expand(trial(n).IDs,length(params.flickerSet),1);
    flip.rsvp(ind) = Expand(trial(n).rsvp,length(params.flickerSet),1);
    for m = 1:2
        if trial(n).targ(m)==1
            targ1 = ind(1)+((trial(n).targInd(m)-1)*length(params.flickerSet));
            flip.targs(m,targ1:targ1+length(params.flickerSet)-1) = ones(1,length(params.flickerSet)); end
    end
end

%%%%%%%%%%%%%%%%%%%%
% open window     %
%%%%%%%%%%%%%%%%%%%%

HideCursor;
Priority(9);

%%%% open screen
sampleFace = face{1};
params.backgroundColor = sampleFace(1,1); % color
screen=max(Screen('Screens'));
[win, rect]=Screen('OpenWindow',screen,params.backgroundColor);
if atScanner Screen(win, 'TextSize', 24);
else Screen(win, 'TextSize', 18); end

Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%%%% timing optimization
scan.flipInt = Screen('GetFlipInterval',win);
slack = scan.flipInt/2;
scan.time = params.countDown+(numTrials*length(params.trialSet)*params.timeBase);   % time in seconds
scan.allFlips = (0:params.flickerBase:scan.time);

%%%% scale the stims for the screen
scan.ppd = pi* rect(3) / (atan(scan.screenWidth/scan.viewingDist/2)) / 360;
scan.faceSize = round(params.faceSizeDeg*scan.ppd);                 % in degrees, the size of our faces
scan.fixRad = round(params.fixRadDeg*scan.ppd);
scan.littleFix = round(scan.fixRad*.25);
%scan.fixImRad = round(params.fixImDeg*scan.ppd);

%%% move xc and yc if needed
xc = rect(3)/2+scan.offset(1); % rect and center, with the flixibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2+scan.offset(2);

%%% define grid of stimulus locations
gr = scan.ppd*params.gridSpaceDeg*([1:params.grid]-ceil(params.grid/2));
[X,Y]=(meshgrid(gr));X=X';Y=Y';
%%% sample centers across, then down
scan.centers = [X(:)+xc Y(:)+yc];

%%%%% flipping for the 16ch headcoil
if flipScreen == 1
    scan.centers = repmat([rect(3) rect(4)],length(scan.centers),1)-scan.centers;
    for n = 1:length(face)
        face{n} = flipud(face{n});
    end
end

scan.rects = CenterRectOnPoint([0 0 scan.faceSize scan.faceSize],scan.centers(:,1), scan.centers(:,2));

if flipScreen && ~onLaptop yc = rect(4)-yc; xc = rect(3)-xc; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       experiment                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% start recording the response
scan.pressFlip = [];

%%% initial window - wait for backtick
Screen(win, 'DrawText', 'Waiting for G Keypress.', 10,10,params.textColor);
DrawFormattedText(win,params.task, 'center', 'center', params.taskColor,[],flipScreen,flipScreen);
Screen(win, 'Flip', 0);

%%%% for CNI
if scan.trigger > 0
    while 1
        KbTriggerWait(KbName('g'), laptop);
        [status, ~] = startScan;
        if status == 0
            break
        else
            DrawFormattedText(win, 'Trigger Failed!', 'center', 'center', params.textColor,[],flipScreen,flipScreen);
            Screen('Flip', win);
        end
    end
else KbTriggerWait(KbName('g'), laptop);
end

%%%% in both CNI/VUIIS cases, this will record keypresses
responseKeys = zeros(1,256);
responseKeys(KbName({'1!';'2@';'3#';'4$';'5%';'6^'}))=1;
KbQueueCreate(scanner,responseKeys);
KbQueueStart();

%%%%%%% START task TASK/FLIPPING
for n = 1:(length(scan.allFlips)-1)
    tr = flip.trials(n);
    
    if flip.IDs(n) > 0 && flip.flicker(n)>0 % draw a face
        if condition(trial(tr).cond).stim == 1 % inv
            f = flipud(face{flip.IDs(n)});
        else f = face{flip.IDs(n)}; end
        fTex = Screen('MakeTexture',win,f);
        
        Screen('DrawTexture', win, fTex,[],scan.rects(condition(trial(tr).cond).pos,:));
        Screen('Close', fTex);
    end
    
    % constant fixation point
    Screen('FillOval', win,[180 0 0], [xc-scan.fixRad yc-scan.fixRad xc+scan.fixRad yc+scan.fixRad]); % outer fixation ring
    Screen('FillOval', win,[255 255 255], [xc-scan.fixRad+params.cuePix yc-scan.fixRad+params.cuePix xc+scan.fixRad-params.cuePix yc+scan.fixRad-params.cuePix]); % white fixation disk
    
    if tr<0 % countdown
        [w,h] = RectSize(Screen('TextBounds',win,num2str(-1*tr)));
        if onLaptop h = -h; end
        DrawFormattedText(win,num2str(-1*tr),xc-w/2,yc-h/2, [0 0 200],[],flipScreen,flipScreen);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% FOR DEBUGGING, MARK ALL TARGETS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     if  flip.targs(1,n) == 1 % letter target
        %         Screen('FillOval', win,[0 0 255 255*.3], [xc-scan.fixRad-5 yc-scan.fixRad-5 xc+scan.fixRad+5 yc+scan.fixRad+5]); % black fixation ring
        %     end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    elseif flip.rsvp(n) > 0  && flip.flicker(n)>0 % draw the rsvp task
        
        % RSVP letter/digit
        [w,h] = RectSize(Screen('TextBounds',win,params.str(flip.rsvp(n))));
        if onLaptop h = -h; end
        DrawFormattedText(win,params.str(flip.rsvp(n)), xc-w/2, yc-h/2,params.taskColor,[],flipScreen,flipScreen);
    elseif flip.rsvp(n) == 0
        Screen('FillOval', win,[180 0 0 255*.5], [xc-3 yc-3 xc+3 yc+3]); % small fixation dot
        Screen('FillOval', win,[255 255 255], [xc-1 yc-1 xc+1 yc+1]); % small fixation dot
    end
    
    %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if n == 1 [VBLT scan.startRun FlipT missed] = Screen(win, 'Flip', 0);
        scan.flipTime(n) = scan.startRun;
    else [VBLT scan.flipTime(n) FlipT missed] = Screen(win, 'Flip', scan.startRun + scan.allFlips(n) - slack);end
    
    % listen for response  - correct if you respond to previous 3 letters
    [pressed, firstPress]= KbQueueCheck();
    if pressed == 1 scan.pressFlip = [scan.pressFlip n];
        KbQueueFlush();
    end
    
end
%%%% to show the very last flip screen for its 500ms
[VBLT scan.flipTime(n+1) FlipT missed] = Screen(win, 'Flip', scan.startRun + scan.allFlips(length(scan.allFlips)) - slack);

%%%%%%%%%%%%%%%%%%
% done! wrap up  %
%%%%%%%%%%%%%%%%%%

scan.runTime = GetSecs - scan.startRun;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% performance
perf = struct('task',params.task,'hitTr',[],'falseFlip',[]);

% for now, this only deals with hits, no false alarms...
for p = scan.pressFlip
    targetRange = p-params.respWindow:p-1;
    if sum(flip.targs(1,targetRange)>0) % we only care about task 1 in this version
        perf.hitTr = [perf.hitTr flip.trials(targetRange(find(flip.targs(1,targetRange))))]; % this is a hit
    else
        perf.falseFlip = [perf.falseFlip p];
    end
end


trialTargs = reshape([trial.targ],2,length(trial));
% we only care about task #1 in this version
perf.targTrials = find(trialTargs(1,:));
perf.hitTr = unique(perf.hitTr);
perf.hitRate = length(perf.hitTr)/length(perf.targTrials);
perf.missTr = setdiff(perf.targTrials,perf.hitTr);
% perf(t).falseRate = length(perf(t).falseTr)/length(perf(t).targTrials);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval(['save data/fixPRF5_' scan.subj '_run' num2str(scan.runNum) '_' scan.date '.mat params scan trial perf flip condition face']);

KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');

fprintf('%s: hit rate %.2f%%\n',params.task,100*perf.hitRate);
fprintf('False alarm presses: %i\n',length(scan.pressFlip)-length([perf.hitTr]));

fclose all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% quick behavior plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % structure with organization: stim(1:2).pos(1:15).perf
% c = struct('perf',[],'overallPerf',[]);
% [c(1:25).perf] = deal([]);
% stim = struct('cond',{'Inverted' 'Upright' 'Blank'},'pos',c);
%
%
% if ~isempty([perf.hitTr]) % if performance is 0, we assume that this run is just a scanner error, and don't record any part of it
%     hitsMisses = {[trial([perf.hitTr]).cond];[trial([perf.missTr]).cond]};
%     pf = [1 0];
%     for m = 1:length(hitsMisses) %hits, misses
%         for t = hitsMisses{m}
%             if t ==0 stimN = 3; posN = 13; % blanks
%             else stimN = condition(t).stim; posN = condition(t).pos; end % for completeness, plot blank performance at the center
%             stim(stimN).pos(posN).perf = ...
%                 [stim(stimN).pos(posN).perf pf(m)];
%         end
%     end
% end
%
% % messy but meh
% c = 1;
% for n = 1:length(stim)
%     for p = 1:length(stim(n).pos)
%         perfPlot(c).name = stim(n).cond;
%         if ~isempty(stim(n).pos(p).perf)
%             stim(n).pos(p).overallPerf = mean(stim(n).pos(p).perf);
%         else stim(n).pos(p).overallPerf = NaN; end
%     end
%     perfPlot(c).vect = [stim(n).pos.overallPerf];
%     perfPlot(c).mat = reshape(perfPlot(c).vect,5,5)';
%     c = c+1;
% end
%
% niceFig([.1 .1 .8 .6],18);
% for c = 1:3
%     subplot(1,3,c)
%     plotInSpace(perfPlot(c).mat,'Behavior',[perfPlot(c).name ': ' num2str(nanmean(perfPlot(c).vect)) ' hit rate'],1,[0 1]);
% end