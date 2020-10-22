% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs =prfSubjs;
expt = 'fixPRF';
tests = {'Y' 'X' 'eccen' 'size' 'gain' 'r2'}; % can be parname (Y,X,sd,gain,exp,shift) or pRF.read value (r2,size,eccen,gain)
whichM = 'mean'; % mean or median

r2cutoff = 'r2-20'; %'perc-50';%'r2-50';%   %    %or 'r2-20'    % cutoff for vox selection
whichANOVA = 'face'; % 'EVC' or 'face';

ROIs= standardROIs(whichANOVA);%{'V1' 'hV4' 'IOG_faces' 'pFus_faces' 'mFus_faces'};
fitSuffix = '';

txtName = [whichANOVA '-' r2cutoff];

whichStim = 'outline';
whichModel = 'kayCSS';
hems = {'rh' 'lh'};

factNames = {'hem' 'ROI' 'condition'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hem = struct;

for h = 1:length(hems)
    pF = ['prfSets/fixPRF_kayCSS_outline_' hems{h} '_r2-20.mat']; % pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems);
    load(pF); fprintf('...%s\n',pF);
    hem(h).subj = subj;
end

ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ttest the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parNames = roi(1).fits(1).parNames;

% data formatting for: RMAOV33.m
% each row = [data fact1# fact2# fact3# subj#];
checkDir([pwd '/stats/']);
    fid = fopen([pwd '/stats/ANOVA3_' txtName '_' whichModel '_' whichStim '_' whichM '.txt'],'w+');
fprintf(fid,'\n**************\n%s, %s\n**************\n',whichM, strTogether(ROIs));
fprintf(fid,'pRF file: %s\n',pF);

fprintf('\n**************\n%s, %s\n**************\n',whichM, strTogether(ROIs));


for t = 1:length(tests)
    test = tests{t};
    anovaData = [];
    testNum = cellNum(test,parNames);
    
    rmSubjs = [];
    for h = 1:length(hems)
        for r = 1:length(ROIs)
            sData = [];
            for s = 1:length(subjs)
                for c = 1:length(hem(1).subj(1).roi(1).fits)
                    % grab voxels
                    vox = hem(h).subj(subjNum(s)).roi(ROInum(r)).fits(c).vox;
                    if ~isempty(vox)
                    if ~isempty(testNum)
                        pars = vertcat(vox.params);
                        eval(['sData = nan' whichM '(pars(:,testNum));']);
                    else
                        eval(['sData = nan' whichM '([vox.' test ']);']);
                    end
                    
                    % an = struct('factor',{'hem' 'ROI' 'condition'},'levels',{hems ROIs {'inverted' 'upright'}});
                    anovaData = [anovaData; sData h r c subjNum(s)];
                    %else error('Missing data! Can''t run this ANOVA function!'); end
                    else rmSubjs(end+1) =  subjNum(s); if c == 1 fprintf('Missing data in %s %s-%s!\n', subjs{s}, hems{h}, ROIs{r});
                        end
                    end
                end
            end
        end
    end
    
    % check for missing data and remove those subjects from the comparison
    [anovaData] = anova_rmSubjs(anovaData,rmSubjs);
    % run anova
    result = rmAnova3(anovaData,factNames,0);
    %print output
    anova3_text(fid,result,test)
    
end
fclose(fid);