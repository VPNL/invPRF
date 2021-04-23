% 2-way anova (bilat ROIs)

clear all; close all;

subjs =prfSubjs;
expt = 'fixPRF';
tests = {'Y' 'X' 'eccen' 'size' 'gain' 'r2' 'sd' 'exp'}; % can be parname (Y,X,sd,gain,exp,shift) or pRF.read value (r2,size,eccen,gain)
whichM = 'mean'; % mean or median

r2cutoff = 'r2-20';
whichTest = 'face';

ROIs= standardROIs(whichTest);

txtName = [whichTest '-' r2cutoff];

whichStim = 'outline';
whichModel = 'kayCSS';%'cssExp1';
hems = {'rh' 'lh'};

factNames = {'ROI' 'condition'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pF =  ['prfSets/fixPRF_' whichModel '_outline_' hemText(hems) '_r2-20.mat'];
    load(pF); fprintf('%s\n...',pF);

ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ttest the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parNames = roi(1).fits(1).parNames;

% data formatting for: RMAOV33.m
% each row = [data fact1# fact2# fact3# subj#];
checkDir([dirOf(pwd) '/stats/' expt]);
    fid = fopen([dirOf(pwd) '/stats/' expt '/ANOVA2_' txtName '_' whichModel '_' whichStim '.txt'],'w+');
fprintf(fid,'\n**************\n%s, %s\n**************\n',whichM, strTogether(ROIs));
fprintf(fid,['pRF file: %s\nROIs: ' repmat('%s ',1,length(ROIs)) '\n'],pF,ROIs{:});

fprintf('\n**************\n%s, %s\n**************\n',whichM, strTogether(ROIs));


for t = 1:length(tests)
    test = tests{t};
    anovaData = [];
    testNum = cellNum(test,parNames);
    
    rmSubjs = [];
        for r = 1:length(ROIs)
            sData = [];
            for s = 1:length(subjs)
                for c = 1:length(subj(1).roi(1).fits)
                    % grab voxels
                    vox = subj(subjNum(s)).roi(ROInum(r)).fits(c).vox;
                    if ~isempty(vox)
                   
                     if ~isempty(testNum)
                        pars = vertcat(vox.params);
                        eval(['sData = nan' whichM '(pars(:,testNum));']);
                    else
                        eval(['sData = nan' whichM '([vox.' test ']);']);
                    end
                    
                    % an = struct('factor',{'hem' 'ROI' 'condition'},'levels',{hems ROIs {'inverted' 'upright'}});
                    anovaData = [anovaData; sData r c subjNum(s)];
                    %else error('Missing data! Can''t run this ANOVA function!'); end
                    else rmSubjs(end+1) =  subjNum(s); if c == 1 fprintf('Missing data in %s %s-%s!\n', subjs{s}, ROIs{r});
                        end
                    end
                end
            end
        end

    
    % check for missing data and remove those subjects from the comparison
    
    [anovaData] = anova_rmSubjs(anovaData,rmSubjs);
    % % function stats = rm_anova2(Y,S,F1,F2,FACTNAMES)
    result = rm_anova2(anovaData(:,1),anovaData(:,end),anovaData(:,2),anovaData(:,3),factNames);
    anova2_text(fid,result,test);
    
end
fclose(fid);

if onLaptop playSound; end