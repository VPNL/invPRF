% an offshoot of the setPRFs procedure that quantifies how many voxels are
% trimmed at various stages; does not actually save a prfset
%%% **** this is here for illustrative purposes only! it requires all of
%%% the original fitting files, which are too cumbersome to release at this
%%% point. 

%clear all; %close all;

info.subjs = prfSubjs;
info.task = '';
info.expt = 'fixPRF';
info.setNotes = '.1 min size + 10.9 max size';

info.whichCutoff ='r2';%'perc';%  % 'perc' or 'r2'
info.cutoffVal = 20;      % cutoff for vox selection - either minR2 or percentile
info.whichPerc = 'min';   % min or mean

if ~strcmp(info.whichCutoff,'perc') info.whichPerc = []; end

info.ROIs= standardROIs;

info.whichStim = 'outline';
info.whichModel ='kayCSS';
info.hems = {'lh' 'rh'};
info.fitSuffix = ''; 

saveFits = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pRF loading - collapse across subjects    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear fits; clear roi;

for r = 1:length(info.ROIs)
    init = 1;
    for s = 1:length(info.subjs)
        sinit = 1;
        [session, numRuns] = vpnlSessions(info.expt,info.subjs{s}); % OPTIONAL: SESSNUM, TASK
        
        for h = 1:length(info.hems)
            
            thisROI = [info.hems{h} '_' info.ROIs{r}];
            [~, fitsName] = fitsDirs(dirOf(pwd),info.expt,session,info.whichStim,vpnlROI(thisROI,session(1:2),info.expt),info.whichModel,info.fitSuffix);
            if exist(fitsName)>0
                load(fitsName);
                
                % add hem text
            for c = 1:length(fits)
                for v = 1:length(fits(c).vox)
                        fits(c).vox(v).hem = h;
                        fits(c).vox(v).stim = fits(c).stim; 
                end
            end
                
                if sinit
                    subj(s).roi(r).fits = fits; sinit = 0;% keep subjects in distinct fields
                else
                    for c = 1:length(fits)
                        subj(s).roi(r).fits(c).vox = [subj(s).roi(r).fits(c).vox fits(c).vox];
                    end
                end
                
                if init
                    roi(r).fits = fits; init = 0; % across-subject struct
                else
                    for c = 1:length(fits)
                        roi(r).fits(c).vox = [roi(r).fits(c).vox fits(c).vox];
                    end
                end
               
            else fprintf('Missing fits in: %s...\n',fitsName);
               saveFits = 0;  
            end
        end % hems
    end %subjs
end %rois


%if saveFits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% translate fit params into more meaningful terms for invPRF expt              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = 1:length(info.ROIs)
    % across subjects
    for c = 1:length(roi(r).fits)
        if isfield(fits,'expN')
            roi(r).fits(c).vox = readPRFs(roi(r).fits(c).vox,fits(1).ppd,fits(1).res, fits(1).expN);
        else    roi(r).fits(c).vox = readPRFs(roi(r).fits(c).vox,fits(1).ppd,fits(1).res,1);end
        % each subject
        for s = 1:length(info.subjs)
            if isfield(fits,'expN')
                subj(s).roi(r).fits(c).vox = readPRFs(subj(s).roi(r).fits(c).vox,fits(1).ppd,fits(1).res, fits(1).expN);else
                subj(s).roi(r).fits(c).vox = readPRFs(subj(s).roi(r).fits(c).vox,fits(1).ppd,fits(1).res,1);end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % which voxels are we plotting? trim edge values, R2 cutoff                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % each subject
    for s = 1:length(info.subjs)
        
        % trimming procedure for this subject
        plotVox = 1:length(subj(s).roi(r).fits(1).vox);
        subj(s).roi(r).voxcount = length(plotVox);
        
        for v = 1:length(subj(s).roi(r).fits(1).vox)
        %%% routine voxel trimming first
                if  ~strcmp(info.expt,'nhp')
                    for c = 1:length(subj(s).roi(r).fits) % for nhp data we only do r2 triming, no additional trimPRFs procedure
                    if trimPRFs2(subj(s).roi(r).fits(c).vox(v),fits(1).ppd,fits(1).res)
                        plotVox(find(plotVox==v)) = []; end
                end
                end
        end
        
        subj(s).roi(r).initTrim = subj(s).roi(r).voxcount - length(plotVox);
        for c = 1:length(subj(s).roi(r).fits)
        subj(s).roi(r).fits(c).vox = subj(s).roi(r).fits(c).vox(plotVox);
        end
        %%%% then, thresholding within subject
        
        
        %%%% determine threshold
            % threshold fits by r2
            clear r2s;
                  if strcmp(info.whichCutoff,'r2') 
                        info.cutoff(r,s)=info.cutoffVal;
                        
            % threshold fits by percentile            
                    elseif strcmp(info.whichCutoff,'perc') 
                        % aggregate r2 across conditions
                        for c = 1:length(subj(s).roi(r).fits)
                        r2s(c,:) = [subj(s).roi(r).fits(c).vox.r2];
                        end
                        eval(['r2s = ' info.whichPerc '(r2s);']);
                        info.cutoff(r,s) = prctile(r2s,info.cutoffVal);
                  end
                  
        %%%% apply threshold          
        plotVox = 1:length(subj(s).roi(r).fits(1).vox);
        for v = 1:length(subj(s).roi(r).fits(1).vox)
          for c = 1:length(subj(s).roi(r).fits) % for nhp data we only do r2 triming, no additional trimPRFs procedure
                    if subj(s).roi(r).fits(c).vox(v).r2 < info.cutoff(r,s) 
                        plotVox(find(plotVox==v)) = []; end      
          end
        end %vox
        subj(s).roi(r).r2Trim = length(subj(s).roi(r).fits(1).vox) - length(plotVox);
        
        %%% pick voxels for this subject & add to general struct
        for c = 1:length(subj(s).roi(r).fits)
            info.numVox(r,s) = length(plotVox);
            subj(s).roi(r).fits(c).vox = subj(s).roi(r).fits(c).vox(plotVox);
       
        if s == 1
            roi(r).fits(c).vox = subj(s).roi(r).fits(c).vox;
        else  roi(r).fits(c).vox = [roi(r).fits(c).vox subj(s).roi(r).fits(c).vox];
        end
        end % condition
        
    end % subj
    if strcmp(info.whichCutoff,'perc')
        fprintf(['%s cutoffs: ' repmat('%.2f ',1,length(info.subjs)) '\n'],...
            info.ROIs{r},info.cutoff(r,:)); 
        fprintf(['%s numVox: ' repmat('%.2f ',1,length(info.subjs)) '\n'],...
            info.ROIs{r},info.numVox(r,:)); end
end % roi

%%% aggregate trim info
for r = 1:length(info.ROIs)
    for s = 1:length(subj)
            initvox(r,s) = subj(s).roi(r).voxcount; % init count
            trim(r,s,1) = subj(s).roi(r).initTrim; 
            trim(r,s,2) = subj(s).roi(r).r2Trim;
            trim(r,s,3) = length(subj(s).roi(r).fits(c).vox); % final count
            if ~isequal(initvox(r,s),sum(trim(r,s,:)))
                fprintf('Not adding up! subj %s, %s...\n',info.subjs{s},info.ROIs{r});
            end
    end
end

initvox_grp = sum(initvox,2)
trim_grp = squeeze(sum(trim,2))
%isequal(initvox_grp,sum(trim_grp,2))

prop = trim_grp./initvox_grp

%bilatFile = pRFfile(dirOf(pwd),info.expt,[info.whichCutoff '-' num2str(info.cutoffVal)],info.whichStim,info.whichModel,info.hems,info.fitSuffix,info.task);
%fprintf('Saving %s...\n',bilatFile);
%save(bilatFile,'roi','subj','info');
barh(prop,'stacked')
legend({'init-trim','r2-trim','incl'});yticklabels(info.ROIs);

%save('trim_info.mat','initvox_grp','trim_grp','prop','initvox','trim','info');