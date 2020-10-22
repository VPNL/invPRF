% prepare css data to be fit via sherlock functions (stanford cluster
% compute resources)

clear all; close all;

whichStim = 'outline';

expt = 'fixPRF';
mainDir = dirOf(pwd);

ROIs = standardROIs;%{'hV4'};%('face');%['V2' 'V3' standardROIs('face')];%{'hV4'};%standardROIs;%%standardROIs('face+');%
hems = {'rh' 'lh'};%,
modelSuffix ='kayCSS';%'%cssExpN';% 'inflipCSSn';% 'intempCSSn';%

subjs = prfSubjs;
getBetas = 0;   % if 1 (or we haven't done it before), will re-compute glm betas
whichDT = 4;    % specify dataType number to work with for betas (default = 4 for sliced, mc data)

for s = 1:length(subjs)
    [sessions{s},nScans] = vpnlSessions(expt,subjs{s});
    stimFile = [dirOf(pwd) 'cssFit/stims/' sessions{s} '_condAvg_' whichStim '.mat'];
    if ~exist(stimFile)
        fprintf('Stim file %s is not specified correctly!Attempting now...\n',stimFile);
        try
            fixPRF_writeFiltIms(whichStim,sessions{s},mainDir,[1:nScans]);
        catch
            error(sprintf('Stim file %s could not be created.\n',stimFile)); end
    end
    
    
    cd([mainDir '/cssFit']);
    
    % specify the index of the voxel to fit
    for h = 1:length(hems)
        for r = 1:length(ROIs)
            % define the final model specification.
            hem = hems{h};
            thisROI= vpnlROI([hem '_' ROIs{r}],subjs{s});
            [dataName, fitsName] = fitsDirs(dirOf(pwd),expt,sessions{s},whichStim,thisROI,modelSuffix);
            
            %%%%% get betas from mrVista
            if ~exist(dataName) || getBetas == 1
                fprintf('Computing betas for %s...\n',dataName);
                recomp = 1;
                outputFile = cssPrep_getBetas(expt,thisROI,sessions{s},numScans,whichDT,recomp);
                fprintf('Done with %s betas, saved to %s...\n',thisROI,outputFile);
            else
                fprintf('Won''t recompute betas for %s...\n',dataName);
            end
        end
    end
end

setpref('VISTA', 'recomputeVoxData',1);

%cd([mainDir '/' expt]);

%%%%%% generate a file that will run this fit configuration
configFile = ['sherlock/run_' modelSuffix '_' whichStim '.m'];
checkDir(dirOf(configFile));
fid = fopen(configFile,'w');
fprintf(fid,'addpath(genpath(''~/utils''));\n');

fprintf(fid,'t = getenv(''SLURM_ARRAY_TASK_ID'');\nt = str2num(t)+1;\n');
fprintf(fid,['sessions = {' repmat('''%s'' ',1,length(sessions)) '};\n'],sessions{:});
fprintf(fid,['ROIs = {' repmat('''%s'' ',1,length(ROIs)) '};\n'],ROIs{:});
fprintf(fid,['hems = {' repmat('''%s'' ',1,length(hems)) '};\n'],hems{:});

fprintf(fid,'d(:,1) = repmat([1:length(sessions)],1,length(ROIs)*length(hems));\n');
fprintf(fid,'d(:,2) = Expand(repmat([1:length(ROIs)],1,length(hems)),length(sessions),1);\n');
fprintf(fid,'d(:,3) = Expand([1:length(hems)],length(ROIs)*length(sessions),1);\n');

fprintf(fid,'cssFit_sherlock(''%s'',''%s'',sessions{d(t,1)},''%s'', ROIs{d(t,2)}, hems{d(t,3)});\n',modelSuffix,expt,whichStim);
fclose(fid);

%%%%% generate a file that will bash the fits
bashFile = ['sherlock/submit_' modelSuffix '_' whichStim '.sh'];
fid = fopen(bashFile,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'#SBATCH --job-name=%s\n',modelSuffix);
fprintf(fid,'#SBATCH --array=0-%d\n',length(ROIs)*length(hems)*length(sessions)-1);
fprintf(fid,'#SBATCH --time=24:00:00\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=12\n#SBATCH --mem-per-cpu=2GB\n#SBATCH -p hns,normal\n');
fprintf(fid,'module load matlab/R2017a\n');
fprintf(fid,'matlab -nodisplay < run_%s_%s.m\n',modelSuffix,whichStim);
fclose(fid);

%%%%%% output the commands that will transfer this over to sherlock (to be run from the terminal)
fprintf('\n\n\n****Sherlock Command Line:****\n');
fprintf('scp %scssFit/sherlock/*%s* sonia09@login.sherlock.stanford.edu:~\n',mainDir,modelSuffix); % copy batch & run file to sherlock

%%% this works fairly poorly now that we're fitting all sessions at once...
fprintf('****Sherlock Grab Fits:****\n');
[~, fitsName] = fitsDirs(dirOf(pwd),expt,sessions{s},whichStim,vpnlROI([hems{1} '_' ROIs{1}],'TH'),modelSuffix);
checkDir(dirOf(fitsName));
fprintf('scp -r sonia09@dtn.sherlock.stanford.edu:/scratch/users/sonia09/fits/%s  %s\n',expt,dirOf(fitsName,3)); % copy data from sherlock - all data for this experiment
for s = 1:length(sessions)
    fprintf('scp -r sonia09@dtn.sherlock.stanford.edu:/scratch/users/sonia09/fits/%s/%s  %s\n',expt,sessions{s},dirOf(fitsName,2)); % copy data from sherlock - this subj
end

fprintf('scp -r %ssonia/utils sonia09@login.sherlock.stanford.edu:~\n',dirOf(pwd,1)); % copy utils folder
fprintf('scp -r %ssonia/utils sonia09@dtn.sherlock.stanford.edu:\n',dirOf(pwd,1)); % copy utils folder

%clear all;
