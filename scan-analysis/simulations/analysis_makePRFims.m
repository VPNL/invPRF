% helper function for prfRec simulation that creates array of all pRFs in
% image space; this cuts down on computational costs needed to generate
% these for every simulation. subsequent code will just load & subsample
% these pRFs
% by default, goes through subject-by-subject so that we have

clear all; close all;

expt = 'compPRF'; % assumes all subject in pRFset
minR2 = ['r2-20'];          % cutoff for vox selection
ROIs= standardROIs([1 4 5 6 7]);%{'mFus_faces'};%{'IOG_faces' 'pFus_faces' 'V1' 'pSTS_faces'}; %{'IOG_faces' 'pFus_faces' 'mFus_faces' 'pSTS_faces'}; % if more than one ROI, combine voxel positions for them

whichStim = 'photo';%'outline';%'internal';%
whichModel = 'kayCSS';%'cssShift';%
suffix= '';
hems = {'rh' 'lh'};
recomp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in data                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now we load in the data from both hemispheres, and threshold across
pRFset = pRFfile([raid 'invPRF/'],expt,minR2,whichStim,whichModel,hems,'');
load(pRFset);

vis.ppd = roi(1).fits(1).ppd;
vis.res = roi(1).fits(1).res+vis.ppd; % give a little buffer around our image/make indexing easier

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make coverage map                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[XX,YY]=meshgrid(-vis.res/2:vis.res/2);
checkDir([raid 'invPRF/fixPRF/pRFims/']);

if containsTxt(suffix,'half') mult = .5; else mult = 1; end

for rr = 1:length(ROIs)
    
    fprintf('** Starting %s ...\n',ROIs{rr});
    tic
    r = cellNum(ROIs{rr},info.ROIs);
    for c = 1:length(roi(1).fits)
        imFile = [raid 'invPRF/fixPRF/prfIms/' ROIs{rr} '_' fileName(pRFset) '_cond' num2str(c) suffix '.mat'];
        if ~exist(imFile) || recomp
            cPRFs = []; vis.subjInd = [];
            for s = 1:length(info.subjs)
                sPRFs = [];
                fprintf('Starting %s %s cond %d (%d voxels)...\n',ROIs{rr},info.subjs{s},c,length(subj(s).roi(r).fits(1).vox));
                %tic
                parfor v = 1:length(subj(s).roi(r).fits(c).vox)
                    sPRFs(v,:,:) = PRF(XX,YY,subj(s).roi(r).fits(c).vox(v).XYdeg(1)*vis.ppd,...
                        subj(s).roi(r).fits(c).vox(v).XYdeg(2)*vis.ppd,subj(s).roi(r).fits(c).vox(v).size*vis.ppd*mult);
                end
                %toc
                cPRFs = [cPRFs; sPRFs];
                vis.subjInd = [vis.subjInd repmat(s,1,size(sPRFs,1))];
            end
            %             eval(['cond' num2str(c) '=cPRFs;']);
            fprintf(['Saving ' imFile '...\n']);
            save(imFile,'vis','cPRFs');
            toc;
        else fprintf(['Won''t recompute ' imFile '\n']); end
    end
end