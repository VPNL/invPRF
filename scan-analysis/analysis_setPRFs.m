% for the plotCode update, we generate and save a struct of pRF fits for
% experiments that can be efficiently used and manipulated by other

%clear all; %close all;

info.subjs = prfSubjs;%prfSubjs;%{'george'};%
info.task = '';
info.expt = 'fixPRF';%'fixPRF';%'nhp';%
info.setNotes = '.05 min size + 10.95 max size';

info.whichCutoff ='r2';%'perc';%  % 'perc' or 'r2'
info.cutoffVal = 20;      % cutoff for vox selection - either minR2 or percentile
info.whichPerc = 'min';   % min or mean

if ~strcmp(info.whichCutoff,'perc') info.whichPerc = []; end

info.ROIs= standardROIs(1:7);%{'V1' 'hV4' 'mFus_faces'};%{'PL' 'ML'};%

info.whichStim = 'outline';%'outline';%'internal';%'internal';%'edge';%';% edge'%'
info.whichModel ='kayCSS';%'inflipCSSn';%'intempCSSn';%'cssExpN';%'kayCSS';% .  %'cssShift';%
info.hems = {'lh' 'rh'};% now will automatically also save rh- and lh- versions of the file
info.fitSuffix = ''; %'_orig';%

saveFits = 1;

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


if saveFits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% translate fit params into more meaningful terms for invPRF expt              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = 1:length(info.ROIs)
    % across subjects
    for c = 1:length(roi(r).fits)
        if isfield(fits,'expN')
            roi(r).fits(c).vox = readPRFs(roi(r).fits(c).vox,fits(1).ppd,fits(1).res, fits(1).expN);
        else    roi(r).fits(c).vox = readPRFs(roi(r).fits(c).vox,fits(1).ppd,fits(1).res,[]);end
        % each subject
        for s = 1:length(info.subjs)
            if isfield(fits,'expN')
                subj(s).roi(r).fits(c).vox = readPRFs(subj(s).roi(r).fits(c).vox,fits(1).ppd,fits(1).res, fits(1).expN);else
                subj(s).roi(r).fits(c).vox = readPRFs(subj(s).roi(r).fits(c).vox,fits(1).ppd,fits(1).res,[]);end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % which voxels are we plotting? trim edge values, R2 cutoff                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % each subject
    for s = 1:length(info.subjs)
        
        % trimming procedure for this subject
        plotVox = 1:length(subj(s).roi(r).fits(1).vox);
        for v = 1:length(subj(s).roi(r).fits(1).vox)
        %%% routine voxel trimming first
                if  ~strcmp(info.expt,'nhp')
                    for c = 1:length(subj(s).roi(r).fits) % for nhp data we only do r2 triming, no additional trimPRFs procedure
                    if trimPRFs2(subj(s).roi(r).fits(c).vox(v),fits(1).ppd,fits(1).res)
                        plotVox(find(plotVox==v)) = []; end
                end
                end
        end
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

bilatFile = pRFfile(dirOf(pwd),info.expt,[info.whichCutoff '-' num2str(info.cutoffVal)],info.whichStim,info.whichModel,info.hems,info.fitSuffix,info.task);
fprintf('Saving %s...\n',bilatFile);
save(bilatFile,'roi','subj','info');

%%% after implementing perc, we will automatically save hem-specific files
%%% so that the same voxels are used in bilat and hem-sorted analyses
hems = info.hems;
bilatROI = roi; bilatSubj = subj;
for h = 1:length(hems)
    info.hems = hems(h);
    for r = 1:length(info.ROIs)
        for c = 1:length(roi(r).fits)
        roi(r).fits(c).vox = bilatROI(r).fits(c).vox(find([bilatROI(r).fits(c).vox.hem]==h));
        for s = 1:length(subj)
            subj(s).roi(r).fits(c).vox = bilatSubj(s).roi(r).fits(c).vox(find([bilatSubj(s).roi(r).fits(c).vox.hem]==h));
            info.hemVox(r,s,h) = length(subj(s).roi(r).fits(c).vox);
        end
        end
    end
        
saveFile = pRFfile(dirOf(pwd),info.expt,[info.whichCutoff '-' num2str(info.cutoffVal)],info.whichStim,info.whichModel,info.hems,info.fitSuffix,info.task);
fprintf('Saving %s...\n',saveFile);
save(saveFile,'roi','subj','info');
end

% plots with some information about these voxels
f = niceFig([.3 .3 .6 .6]); 
subplot(3,1,1); title('R2 cutoffs'); bar(info.cutoff,'grouped'); set(gca,'xticklabel',info.ROIs,'box','off'); legend(info.subjs);
subplot(3,1,2:3); title('Voxel Counts'); plotBarStackGroups(info.hemVox,info.ROIs,hems)
superTitle(fileName(bilatFile),.05);
niceSave([dirOf(bilatFile) 'ims/'],fileName(bilatFile),[],[],{'png'}); %close(f);
% 
% figure(modelComp); subplot(2,3,m);
if strcmp(info.whichCutoff,'perc') 
niceFig([.4 .4 .2 .3]);niceBars2(info.cutoff','mean',1,info.ROIs);t = title(fileName(bilatFile)); set(t,'interpreter','none'); end %set(gca,'box','off');
% m = m+1;
%niceSave([dirOf(bilatFile) 'ims/'],['mean_' fileName(bilatFile)],[],[],{'png'});

else
    fprintf('** Won''t save a pRF set...\n'); end