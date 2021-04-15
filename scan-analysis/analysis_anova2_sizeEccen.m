% loads & plots model properties (like size x eccen) along the visual heirarchy for cssFit properties
% 3/11/20: THIS IS THE SIZE ANOVA CODE
% 5/6/20: changed size to == 1 sigma/sqrt(n)

clear all; close all;


expt = 'fixPRF';

minR2 = 20;          % cutoff for vox selection
whichANOVA = 'face';
ROIs= [standardROIs(whichANOVA)];%('face');% 'V3' 'hV4' 'IOG_faces' 'pFus_faces','mFus_faces'
saveFigs = 0;

whichStim = 'outline';
whichModel = 'kayCSS';

hems = {'lh' 'rh'};
tests = {'slope' 'intercept'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveDir = [pwd '/stats/'];
checkDir(saveDir);
fid = fopen([saveDir 'ANOVA2_' hemText(hems) '_sizeEccen_' whichANOVA '-' num2str(minR2) '_' whichModel '_' whichStim '.txt'],'w+');

fontSize = 11; titleSize = 14;
for t = 1:length(tests)
    anovaData= []; rmSubjs = [];
    pf = ['prfSets/fixPRF_kayCSS_outline_' hemText(hems) '_r2-20.mat']; % pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems);
    fprintf(fid,'\n**************\n%s\n**************\n',strTogether(ROIs));
    load(pf);fprintf(fid,'pRF file: %s\n',pf);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure 1: vox scatter + fit lines
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for r = 1:length(ROIs)
        if saveFigs niceFig([.1 .1 .8 .8],fontSize,1); end
        
        ROInum = cellNum(ROIs{r},info.ROIs);
        sFit = struct;
        for s = 1:length(subj)
            fits = subj(s).roi(ROInum).fits;
            
            for c = 1:length(fits)
                if c == 1 mult = .25; else mult = 1; end
                
                
                x = [fits(c).vox.eccen]';
                y = [fits(c).vox.size]'/2;
                X = [x ones(length(fits(c).vox),1)];
                
                if ~isempty(X)
                    [sFit(c).h(s,:),R2(s)] = fitl1line(X,y);
                    if saveFigs subplot(2,length(info.subjs)/2,s); 
                    h1 = plot(x,X*sFit(c).h(s,:)','Color',condColors(s,1)*mult); hold on;
                    hold on; scatter(x,y,5,condColors(s,1)*mult,'filled'); hold on; title(['bilat ' ROIs{r} ' subj ' info.subjs{s}]);
                    set(h1,'LineWidth',1); %set(extLine,'LineStyle',':');
                    alpha(.2);
                    end
                    anovaData = [anovaData; sFit(c).h(s,t) r c s];
                else rmSubjs(end+1) = s;
                end
            end
            
            if saveFigs
            xl = xlim; xlim([0.25 6]);  ylim([0.25 6]);
            set(gca,'TickDir','out');
            xlabel('Eccen (dva)','FontSize',fontSize); ylabel(['Size (1*SD/sqrt(N)) (dva)'],'FontSize',fontSize);
            axis square; end
            
        end
    if saveFigs niceSave(saveDir,['anova2data_' hemText(hems) '-' ROIs{r}]); end    
    end
    
    fprintf('\n**************\n%s\n**************\n', strTogether(ROIs));
    
    % check for missing data and remove those subjects from the comparison
    [anovaData] = anova_rmSubjs(anovaData,rmSubjs);
    % run anova
    factNames = {'ROI' 'condition'};
    result = rm_anova2(anovaData(:,1),anovaData(:,end),anovaData(:,2),anovaData(:,3),factNames);
    %print output
    anova2_text(fid,result,tests{t});
    
    fprintf('\n*****\nt-tests by condition:\n****\n');
    %%% ROI-wise ttesting
    for r = 1:length(ROIs)
        rInd = anovaData(:,2) == r;
        rData = anovaData(anovaData(:,2) == r,:);
        for c = 1:2
            sData(:,c) = rData(rData(:,3)==c,1); end
        [H,P,CI,STATS] = ttest(sData(:,1),sData(:,2));
        diff = mean(sData(:,2)-sData(:,1));
        se = std(sData(:,2)-sData(:,1))/sqrt(length(sData));
        if strcmp(tests{t},'X') || strcmp(tests{t},'Y')
            diff = diff /roi(1).fits(1).ppd;
            se = se /roi(1).fits(1).ppd;end
        p(r) = P;
        stats{r} = STATS;
        if H sig = '***'; else sig = ''; end
        fprintf('%s [%s %s:] %s: t(%d)=%.5f, p=%.5f; diff = %.3f (SE=%.3f)\n',sig,hemText(hems),ROIs{r},tests{t},STATS.df,STATS.tstat,P,diff,se);
        fprintf(fid,'%s [%s %s:]  %s: t(%d)=%.5f, p=%.5f; diff = %.3f (SE=%.3f)\n',sig,hemText(hems),ROIs{r},tests{t},STATS.df,STATS.tstat,P,diff,se);
    end
end
playSound;