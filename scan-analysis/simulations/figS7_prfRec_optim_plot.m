% simulates overlap of pRF coverage and mean face features at different
% locations
% annoyingly, this does require PTB...
% _parallel gives us fewer checkpoints/output but should be faster

%clear all; close all;


expt = 'fixPRF'; % assumes all subject in pRFset
minR2 = ['r2-20'];          % cutoff for vox selection
ROIs = {'mFus_faces'};%{'IOG_faces' 'pFus_faces' 'V1' 'pSTS_faces'};%{'IOG_faces' 'pFus_faces' 'mFus_faces' 'pSTS_faces'}; % if more than one ROI, combine voxel positions for them

whichStim = 'outline';%'photo';%'internal';%
whichLoc = 'outline';
whichModel = 'kayCSS';%'cssShift';%
hems = {'rh' 'lh'};

sim.suffix = ''; % if blank, uses original convention of setting prf sd to equal size as 2 sigma/sqrt(n);
sim.numSims = 1000;%50;
sim.drawVox = .8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in data, sort plotVox           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
niceFig([.1 .1 .9 .9]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load make maps of group & indiv optims  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for c = 1:2
    
    for r = 1:length(ROIs)
        if length(ROIs)>1
        subplot(2,length(ROIs),(c-1)*length(ROIs)+r);
        else subplot(1,2,c); end
        
        simName = [raid 'invPRF/fixPRF/behavSim/' whichStim '_' ROIs{r} '_' minR2 '_' num2str(sim.drawVox) 'vox' sim.suffix];
        if ~strcmp(whichLoc,'internal') simName = [simName '_' whichLoc]; end
         
        if ~exist([simName '_cond' num2str(c) '.mat'])
            error(sprintf(['Missing ' [simName '_cond' num2str(c)] '.mat! Run prfRec_optim_monteCarlo.m to compute...\n']));
        else
            load([simName '_cond' num2str(c) '.mat']);
        end
        XY = vertcat(sim.subj([1 2 4 6 10 11
            3 5 7 8 9 12]).bestPos);  % prfSubjs aren't actually in chrono order, but this should fix it
        
        
        %%% plot individual subject positions
        
        prfAxes(sim.ppd,sim.res);
        
        subjCols = [colorinterpolate([condColors(1)*1.5; condColors(2)],6,0); colorinterpolate([condColors(1,1); condColors(2,1)],6,0)];
        %subjCols = subjCols(,:); % prfSubjs aren't actually in chrono order, but this should fix it

        
        for x = 1:size(XY,1)
            s(x) = scatter(XY(x,1),XY(x,2),20,subjCols(x,:),'filled'); hold on; end
        legend(s,prfSubjs([1 2 4 6 10 11 ...
            3 5 7 8 9 12]));
        grid on;
        %%% underlay face image at ideal location
        bestPos = mean(XY); sem = se(XY);
        %hold on; scatter(bestPos(1),bestPos(2),50,[.5 .5 .5],'*'); hold on;
        faceIm = zeros(sim.res+1,sim.res+1); co = CenterRectOnPoint([1 1 (sim.faceSize*sim.ppd) (sim.faceSize*sim.ppd)],...
            round(bestPos(2)*sim.ppd)+sim.res/2+1,round(bestPos(1)*sim.ppd)+sim.res/2+1);
        faceIm(co(1):co(3),co(2):co(4)) = face;
        
        h = imagesc(xlim,ylim,imcomplement(faceIm)); colormap gray;
        uistack(h,'bottom');
        title({ROIs{r};sprintf('Best Location: [%.2f (SE=%.2f) %.2f (SE=%.2f)] deg',bestPos(1),sem(1),bestPos(2),sem(2))});   
    end
end
%superTitle(sprintf('Simulation:,12,.05);
niceSave([dirOf(pwd) 'figures/' expt '/optimPosition/'],['fullSim_' whichStim '_' whichLoc],ROIs,[],{'png' 'svg'}); % just save pngs, since these can be generated pretty quickly
    