clear all; close all;

% experiment and session

plotPars = {'gain' 'r2' 'Y' 'X' 'size' };%{'Y'};%
parTitles = { 'Gain Estim' 'Estimated R^{2}' 'Y Estim.' 'X Estim.' 'Size [2xSD/sqrt(N)] (dva)'};%{'Y estim'};%
whichPlots = [0 1]; % 1 = single random voxel, 2 = summary of all plots
saveFig = 1;

sim.expt = 'fixPRF';
sim.whichModel = 'kayCSS'; sim.whichStim = 'outline'; sim.minR2 = 20;
sim.simSuffix = 'scaleIterNoise';

sim.ROI = 'mFus_faces';
sim.hem = 'rh';

load(pRFfile(dirOf(pwd),sim.expt,sim.minR2,sim.whichStim,sim.whichModel,{sim.hem}));
fits = roi(cellNum(sim.ROI,info.ROIs)).fits; % since we're looking at one ROI at a time here, simplify


load(['iterNoiseSim/' sim.simSuffix '_' sim.expt '_' sim.whichModel '_' sim.hem sim.ROI '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% single vox plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if whichPlots(1) == 1
    v = randsample(1:sim.numVox,1);
    
    %%%%%% plot and save each of the modelfits
    if onLaptop figSize = [.2 .2 .5 .5];
    else figSize = [.2 .2 .3 .4];
    end
    f = niceFig(figSize,12);
    numPlots = [1 4];
    
    %%%%%%%%% plot original betas and modelfit
    subplot(numPlots(1),numPlots(2),1:2);
    
    hbar(1) = niceBars(fits(sim.baseCond).condNums,[sim.base(v).betas],[sim.base(v).sems],condColors(sim.baseCond));
    plot(fits(sim.baseCond).condNums,sim.base(v).modelfit,'b-','LineWidth',2);
    
    hbar(2) = niceBars(fits(sim.compCond).condNums,[sim.comp(v).betas],[sim.comp(v).sems],condColors(sim.compCond));
    plot(fits(sim.compCond).condNums,sim.comp(v).modelfit,'b-','LineWidth',2);
    
    legend([hbar(1:2)],{fits(sim.baseCond).cond fits(sim.compCond).cond},'location','SouthEast','FontSize',10,'box','off');
    
    xlabel('Stimulus Number'); xlim([0 fits(end).condNums(end)+1]);
    ylabel('Simulated Betas');
    
    titleTx{1} = ['Vox ' num2str(v) ', R^{2}=' num2str(round(fits(1).vox(sim.voxInd(v)).r2)) ' ' fits(1).cond ];
    for t = 2:length(fits)
        titleTx{t} = ['R^{2}=' num2str(round(fits(t).vox(sim.voxInd(v)).r2)) ' ' fits(t).cond];
    end
    title(titleTx);
    yl=ylim;
    set(gca,'box','off','FontName','Arial','FontWeight','normal','TickDir','out');
    
    %%%%%%%%% plot new betas and modelfit
    subplot(numPlots(1),numPlots(2),3);
    hbar(1) = niceBars([1:25],[sim.vox(v).simData],zeros(1,25),condColors(3)); hold on;
    hbar(2) = niceBars([1:25],[sim.vox(v).addNoise],zeros(1,25),condColors(4)); hold on;
    plot([1:25],sim.vox(v).modelfit,'b-','LineWidth',2);
    
    legend(hbar(1:2),{['Simulated from ' fits(sim.baseCond).cond] 'Added Noise'},'location','SouthEast','FontSize',10,'box','off');
    ylim(yl); %match to prev plot
    
    xlabel(['Noise Level: ' num2str(sim.vox(v).noiseLevel)]); xlim([0 26]);
    ylabel('Simulated Betas');
    
    title([sim.simSuffix ' Simulation, R^{2}=' num2str(round(sim.vox(v).r2))]);
    
    set(gca,'box','off','FontName','Arial','FontWeight','normal','TickDir','out');
    
    %%%%%% fits in space
    subplot(numPlots(1),numPlots(2),4);
    
    for c = 1:length(fits)
        % PRF size is defined as S/sqrt(N).
        h(c) = plotCoverage(fits(c).vox(sim.voxInd(v)),condColors(c),[],fits(1).ppd,fits(1).res);
    end
    
    h(c+1) = plotCoverage(sim.vox(v),condColors(3),[],fits(1).ppd,fits(1).res);
    for tt = 1:length(fits)
        parText{tt} = [ fits(tt).cond ' Params: [' num2str(round(fits(tt).vox(sim.voxInd(v)).params)) ']'];
    end
    parText{tt+1} = ['Sim Params: [' num2str(round(sim.vox(v).params)) ']'];
    t = title(parText);
    set(t,'visible','on');
    l=legend([h(1:length(fits)+1)],{fits.cond 'Simulation'}); set(l,'box','off','location','Best');
    
    titleText = ['Iter ' num2str(i) ' ' sim.ROI ' Voxel #' num2str(sim.voxInd(v)) ', ' sim.whichStim ' stim, ' sim.whichModel ' model'];
    
    titleText = [titleText ', Session = ' sim.vox(v).session];
    superTitle(titleText,12, .01)
    
    if saveFig == 1
        txt = [sim.hem '_' sim.ROI '_vox' num2str(v) ];
        niceSave([dirOf(pwd) 'figures/' expt '/noiseIter/'],txt,[],[],{'svg'}); % just save pngs, since these can be generated pretty quickly
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% summary plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if whichPlots(2) == 1
    if saveFig && onLaptop figSize = [.1 .1 .9 .9]; else figSize = [.3 .3 .7 .7]; end
    niceFig(figSize);
    numPlots=[1,length(plotPars)];
    
    % remove voxels that weren't simulated correctly
   % failed = find(sim.vox.
    
    for p = 1:length(plotPars)
        parNum = cellNum(plotPars{p},fits(1).parNames);
        
        % simulation, base condition, comparison condition
        sData = vertcat(sim.vox.params);
        bData = vertcat(sim.base.params);
        cData = vertcat(sim.comp.params);
        
        
        if isempty(parNum)
            switch plotPars{p}
                case 'size'
                    bData = (2*(bData(:,3)'./fits(1).ppd))./sqrt(bData(:,5)');
                    cData = (2*(cData(:,3)'./fits(1).ppd))./sqrt(cData(:,5)');
                    sData = (2*(sData(:,3)'./fits(1).ppd))./sqrt(sData(:,5)');
                case 'eccen'
                    bData = sqrt([(bData(:,2)-im.size/2)/fits(1).ppd].^2 + [-(bData(:,1)-im.size/2)/fits(1).ppd].^2);
                    cData = sqrt([(cData(:,2)-im.size/2)/fits(1).ppd].^2 + [-(cData(:,1)-im.size/2)/fits(1).ppd].^2);
                    sData = sqrt([(sData(:,2)-im.size/2)/fits(1).ppd].^2 + [-(sData(:,1)-im.size/2)/fits(1).ppd].^2);
                case 'r2'
                    bData = [sim.base.r2];
                    cData = [sim.comp.r2];
                    sData = [sim.vox.r2];
            end
        else
            sData = sData(:,parNum)';
            bData = bData(:,parNum)';
            cData = cData(:,parNum)'; end
        
        subplot(numPlots(1),numPlots(2),p)
        
        s(1) = scatterCent(bData,sData,condColors(sim.baseCond),...
            [],[],[],12,1,1);
        hold on;
        s(2) = scatterCent(bData,cData,condColors(sim.compCond),...
            fits(sim.baseCond).cond,'Simulation',parTitles{p},12,1,1);
        legend([s(1) s(2)],{'Simulation' fits(sim.compCond).cond})
        
        %title(ROIs{r},'fontSize',titleSize,'interpreter','none','FontWeight','bold');
        
        %subplot(numPlots(1),numPlots(2),p+numPlots(2))
        
        
    end
        superTitle([sim.hem '_' sim.ROI num2str(sim.numVox) 'voxels rescaled + noise iterated'], .97);
        if saveFig == 1
            txt = [sim.hem '_' sim.ROI];
            
            niceSave([dirOf(pwd) 'figures/fixPRF/noiseIter/'],txt,[],[],{'svg' 'png'}); % just save pngs, since these can be generated pretty quickly
        end
end