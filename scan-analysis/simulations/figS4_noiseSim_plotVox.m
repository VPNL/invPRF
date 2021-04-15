clear all; close all;

% experiment and session
saveFig = 1;

sim.expt = 'fixPRF';
sim.whichModel = 'kayCSS'; sim.whichStim = 'outline'; sim.minR2 = 20;
sim.simSuffix = 'scaleOnly';

sim.ROI = 'mFus_faces';
sim.hem = 'rh';

load(pRFfile(dirOf(pwd),sim.expt,sim.minR2,sim.whichStim,sim.whichModel,{sim.hem}));
fits = roi(cellNum(sim.ROI,info.ROIs)).fits; % since we're looking at one ROI at a time here, simplify

simFile = ['iterNoiseSim/' sim.simSuffix '_' sim.expt '_' sim.whichModel '_' sim.hem sim.ROI '_' sim.simSuffix '.mat'];
load(simFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% single vox plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
voxInd = find([sim.vox.iter]);
for v = voxInd
    
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
    
    titleTx{1} = ['Vox ' num2str(v) ', R^{2}=' num2str(round(sim.base(v).r2)) ' ' fits(sim.baseCond).cond ];
    titleTx{2} = ['R^{2}=' num2str(round(sim.comp(v).r2)) ' ' fits(sim.compCond).cond ];
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
    
     % PRF size is defined as S/sqrt(N).
        h(1) = plotCoverage(sim.base(v),condColors(sim.baseCond),[],fits(1).ppd,fits(1).res);
        h(2) = plotCoverage(sim.comp(v),condColors(sim.compCond),[],fits(1).ppd,fits(1).res);
        h(3) = plotCoverage(sim.vox(v),condColors(3),[],fits(1).ppd,fits(1).res);
        
    parText{1} = [fits(sim.baseCond).cond ' Params: [' num2str(round(sim.base(v).params)) ']'];
    parText{2} = [fits(sim.compCond).cond ' Params: [' num2str(round(sim.comp(v).params)) ']'];
    parText{3} = ['Sim Params: [' num2str(round(sim.vox(v).params)) ']'];
    
    t = title(parText);
    set(t,'visible','on');
    l=legend([h(1:3)],{fits(sim.baseCond).cond fits(sim.compCond).cond 'Simulation'}); set(l,'box','off','location','Best');
    
    titleText = [sim.ROI ' Voxel #' num2str(v) ', ' sim.simSuffix];
    
    titleText = [titleText ', Session = ' sim.vox(v).session];
    superTitle(titleText,12, .01)
    
    if saveFig == 1
        txt = [sim.hem '_' sim.ROI '_vox' num2str(v)];
        niceSave([dirOf(pwd) 'figures/fixPRF/noiseIter/voxPlots-' sim.simSuffix '/'],txt,[],[],{'png'}); % just save pngs, since these can be generated pretty quickly
    end
    close(f);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% summary plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

niceFig([.1 .1 .8 .8]);
plotPars = {'X' 'Y' 'size' 'gain' 'r2'};
    numPlots=[2,length(plotPars)];
    
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
        
        s(1) = scatterCent(bData,sData,condColors(3),...
            [],[],[],12,1,1);
        hold on;
        s(2) = scatterCent(bData,cData,condColors(sim.compCond),...
            fits(sim.baseCond).cond,'Comparison',plotPars{p},12,1,1);
        legend([s(1) s(2)],{'Simulation' fits(sim.compCond).cond})
        
        %title(ROIs{r},'fontSize',titleSize,'interpreter','none','FontWeight','bold');
        
        subplot(numPlots(1),numPlots(2),p+numPlots(2))
        
        s(1) = scatterCent(cData,sData,condColors(3),...
            [],[],[],12,1,1);
        hold on;
        s(2) = scatterCent(cData,bData,condColors(sim.baseCond),...
            fits(sim.compCond).cond,'Comparison',plotPars{p},12,1,1);
        legend([s(1) s(2)],{'Simulation' fits(sim.baseCond).cond})
        
        
    end
        superTitle([sim.hem '_' sim.ROI num2str(sim.numVox) ' ' sim.simSuffix], .97);
         if saveFig == 1
        txt = ['summary_' sim.simSuffix '_' sim.hem '_' sim.ROI];
        niceSave([dirOf(pwd) 'figures/fixPRF/noiseIter/voxPlots-' sim.simSuffix '/'],txt,[],[],{'png' 'svg'}); % just save pngs, since these can be generated pretty quickly
    end
if onLaptop playSound; end