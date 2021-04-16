% analysis script for second preprocessing step for fixPRF scan files -
% after they are trimmed and spike remoed, perform detrending

%%% ET NAME: always in the format subj/fixPRF_runNum - will be parsed later
clear all; close all;

subjs = prfSubjs;
saveFig = 1;
polydegs = [0 1 2];          % polynomial degrees to use to detrend

for s = 1%:length(subjs)
    
    for r = 1:10
        % input
        etName = sprintf('%s/fixPRF_%i',subjs{s},r);
        
        %dataFile = [exptDir 'data/' dataName  '.mat'];
        matDir = [pwd '/mats/'];%[raid 'invPRF/eyetracking/mats/'];
        etMat = [matDir etName '_preprocessed.mat'];
        
        detrendFile = [matDir etName '_preprocessed' num2str(length(polydegs)) '.mat'];
        figDir =  [pwd '/figures/detrending/'];
        
        if exist(etMat)
            fprintf('Starting %s...\n',etMat);
            load(etMat);
            
            niceFig([.1 .1 .9 .7]);plotText = {'X position ' 'Y position '};
            
            for p = 1:2
                % make plot
                set(gca,'YDir','reverse'); subplot(2,2,p);
                
                plot(samples(:,1),samples(:,p+1),'LineWidth',2,'Color',condColors(p));
                xlim([0 samples(end,1)]);
                title(['Preprocessed Data - ' plotText{p}]); xlabel('Trial Time (s)'); ylabel(plotText{p});
            end
            
            
            
            % detrend each run
            
            eyex = samples(:,2); eyey = samples(:,3);
            badix = isnan(eyex);
            
            
            % detrend each run
            polymatrix = constructpolynomialmatrix(length(eyex),polydegs); assert(length(eyex)==length(eyey));
            
            X = polymatrix(~badix,:);
            h = inv(X'*X)*X'*eyex(~badix);
            eyex(~badix) = eyex(~badix) - X*h;
            h = inv(X'*X)*X'*eyey(~badix);
            eyey(~badix) = eyey(~badix) - X*h;
            
            % figure;
            % subplot(2,2,1); plot(samples(:,2));
            % subplot(2,2,2); plot(samples(:,3));
            % subplot(2,2,3); plot(eyex);
            % subplot(2,2,4); plot(eyey);
            % superTitle(num2str(r));
            
            samples(:,2) = eyex; samples(:,3) =eyey;
            for p = 1:2
                % make plot
                set(gca,'YDir','reverse'); subplot(2,2,p+2);
                
                plot(samples(:,1),samples(:,p+1),'LineWidth',2,'Color',condColors(p));
                xlim([0 samples(end,1)]);
                title(['Detrended Data - '  plotText{p}]); xlabel('Trial Time (s)'); ylabel(plotText{p});
            end
            
            save(detrendFile,'samples','rate','sampleFile','startTrim','endTrim','polydegs');
            
            superTitle([etName ' 0-' num2str(length(polydegs)-1) 'deg polynomial detrending'],12,.05);
            if saveFig
                niceSave([figDir subjs{s} '/'],['poly' num2str(length(polydegs)) '_fixPRF_run' num2str(r)]);
                close all;
            end
            
        else
            fprintf('*** Skipping %s...\n',etMat);
        end
    end
    if onLaptop playSound; end
end