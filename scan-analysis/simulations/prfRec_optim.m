% simulates overlap of pRF coverage and mean face features at different
% locations
% streamlined version: swithces to proportion of voxels to estimate,
% eliminates pooled group estimate, eliminates pooled ROI estimate, does
% multiple ROIs in forloop

clear all; close all;

expt = 'fixPRF'; % assumes all subject in pRFset
minR2 = ['r2-20'];          % cutoff for vox selection
ROIs = {'mFus_faces'}; % if more than one ROI, combine voxel positions for them

whichStim ='outline';% 'internal';%'photo';%
whichLoc = 'outline'; % outline or internal
whichModel = 'kayCSS';%'cssShift';%
hems = {'rh' 'lh'};
recomp = 1;

sim.numSims = 1000;%50;
sim.drawVox = .8; % now a proportion, vs. absolute number
sim.suffix =  ''; % now can call halfSize correction; leave blank for newer pRFsets

for r = 1:length(ROIs)
    sim.name = [raid 'invPRF/fixPRF/behavSim/' whichStim '_' ROIs{r} '_' minR2 '_' num2str(sim.drawVox) 'vox' sim.suffix];
    if ~strcmp(whichLoc,'internal') sim.name = [sim.name '_' whichLoc]; end
    for c = 1:2
        if ~exist([sim.name '_cond' num2str(c) '.mat']) || recomp
            fprintf(['Starting ' [sim.name '_cond' num2str(c)] '.mat: %s...\n'],datestr(now)); tic;
            sim.subj = [];
            
            % in this version, we load in pre-made coverage ims of all pRFs rather
            % than calculating them using the PRF() functions
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % load coverage maps                   %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            try
                imFile = [raid 'invPRF/fixPRF/prfIms/' ROIs{r} '_' fileName(pRFfile('',expt,minR2,whichStim,whichModel,hems)) '_cond' num2str(c) sim.suffix '.mat'];
                load(imFile)
            catch
                error(sprintf('Missing %s! Check filepath or run analysis_makePRFIms.m\n',imFile)); end
            
            sim.ppd = vis.ppd;% now taking directly from the pre-made coverages
            sim.res = vis.res; % give a little buffer around our image/make indexing easier
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % simulation parameters                %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [XX,YY]=meshgrid(-sim.res/2:sim.res/2);
            sim.steps = -3:.1:3; [X,Y]=(meshgrid(sim.steps));
            sim.degPos = [Y(:) X(:)];
            sim.centers = sim.degPos*sim.ppd+sim.res/2+1;
            sim.faceSize = 3.2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % internal features images             %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if r == 1 load([whichLoc 'Avg.mat']);
                faceSize = sim.faceSize*sim.ppd;
                face = imresize(avgFace,[faceSize faceSize]);
                if c == 1 face = flipud(face); end
            end
            
            for n = 1:length(sim.centers)
                    % make current face im
                    fIm = zeros(sim.res+1,sim.res+1); co = CenterRectOnPoint([1 1 faceSize faceSize],sim.centers(n,1),sim.centers(n,2));
                    fIm(co(1):co(3),co(2):co(4)) = face;
                    faceIm{n} = fIm;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % run indiv simulation                 %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            imshow(face);drawnow;
            tic
            for n = 1:length(unique(vis.subjInd))
                clear draw; draw(sim.numSims) = struct;
                vox = cPRFs(vis.subjInd==n,:,:);
                % tic
                drawVox = round(size(vox,1)*sim.drawVox);
                if size(vox,1)<sim.drawVox drawVox = size(vox,1); end
                for s = 1:sim.numSims
                    covIm = vox(datasample(1:size(vox,1),drawVox,'Replace',true),:,:);
                    
                    % to deal with the issue of plot/imagesc using different coordinate systems
                    covIm = flipud(squeeze(mean(covIm)));
                    
                    for m = 1:length(sim.centers)
                        featCov = faceIm{m}.*covIm;
                        draw(s).result(m) = sum(featCov(:));
                    end
                    
                    draw(s).best = find(draw(s).result==max(draw(s).result));
                    draw(s).bestDegX = [sim.degPos(draw(s).best,2)];
                    draw(s).bestDegY = [sim.degPos(draw(s).best,1)];
                    
                    if mod(s,100)==0
                        fprintf('Subj %d, done with %i draws...\n',n,s);
                        toc
                    end
                end
                %%%% aggregate subj results
                sim.subj(n).bestPos = [mean([draw.bestDegX]) mean([draw.bestDegY])];
                sim.subj(n).sd = [std([draw.bestDegX]) std([draw.bestDegY])];
                save([sim.name '_cond' num2str(c) '.mat'],'draw','sim','face');
                fprintf(['Saving Subj ' num2str(n) ' sim: ' [sim.name '_cond' num2str(c)] '\n']); toc;
            end
        else fprintf(['Wont recompute ' [sim.name '_cond' num2str(c)] '...\n']); end
    end
end