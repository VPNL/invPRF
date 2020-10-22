% simulates overlap of pRF coverage and mean face features at different
% locations

clear all; close all;


expt = 'fixPRF'; 

saveFig = 1;

minR2 = 'r2-20';%perc-50';          % cutoff for vox selection
ROI= 'pFus_faces';

whichStim = 'outline';%'photo';%'internal';%
whichModel = 'kayCSS';%'cssShift';%
hems = {'rh' 'lh'};
cond = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in data, sort plotVox           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;

% now we load in the data from both hemispheres, and threshold across
load(pRFfile([raid 'invPRF/'],expt,minR2,whichStim,whichModel,hems,''));

res = roi(1).fits(1).res; ppd = roi(1).fits(1).ppd;

res = res+ppd*5; % give a little buffer around our image/make indexing easier

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make coverage map                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:length(subj)
    tic
subjConv(s).name = info.subjs{s};

fits = subj(s).roi(cellNum(ROI,info.ROIs)).fits(cond); % one ROI at a time, for now
[X,Y]=meshgrid(-res/2:res/2);

for v = 1:length(fits.vox)
    subjCov(s).cov(v,:,:) = PRF(X,Y,fits.vox(v).XYdeg(1)*ppd,fits.vox(v).XYdeg(2)*ppd,fits.vox(v).size*ppd);
end

% to deal with the issue of plot/imagesc using different coordinate systems
subjCov(s).covIm = flipud(squeeze(mean(subjCov(s).cov))); 
niceFig([.1 .1 .9 .9],12); title(['pRF Coverage x ' fits.cond ': ' info.subjs{s} ', ' num2str(length(fits.vox)) 'vox, ' ROI]);
subplot(1,4,1); imshow(subjCov(s).covIm); colormap(mrvColorMaps('hot'));title([ROI ' coverage']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation parameters                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sim.steps = -3:.1:3; [X,Y]=(meshgrid(sim.steps));
sim.degPos = [Y(:) X(:)];
sim.centers = sim.degPos*ppd+res/2+1; 
sim.faceSize = 3.2;

% load & resize face
load('stims/internalAvg.mat'); if cond == 1 avgFace = flipud(avgFace); end
faceSize = sim.faceSize*ppd;
face = imresize(avgFace,[faceSize faceSize]);

for n = 1:length(sim.centers)
% make current face im
faceIm = zeros(res+1,res+1); co = CenterRectOnPoint([1 1 faceSize faceSize],sim.centers(n,1),sim.centers(n,2));
faceIm(co(1):co(3),co(2):co(4)) = face;

 subplot(1,4,2); imshow(faceIm);

featCov = faceIm.*subjCov(s).covIm; subplot(1,3,3); imshow(featCov);
result(n) = sum(featCov(:));
subjCov(s).result(n) = result(n);
end

% plot results

subjCov(s).best = find(subjCov(s).result==max(subjCov(s).result));
subjCov(s).bestDegX = [sim.degPos(subjCov(s).best,2)];
subjCov(s).bestDegY = [sim.degPos(subjCov(s).best,1)];
fprintf('Subject %s Best Location: [%.2f %.2f] deg\n',subjConv(s).name,subjCov(s).bestDegX,subjCov(s).bestDegY);

faceIm = zeros(res+1,res+1); co = CenterRectOnPoint([1 1 faceSize faceSize],round(subjCov(s).bestDegY *ppd)+res/2+1,round(subjCov(s).bestDegX *ppd)+res/2+1);
faceIm(co(1):co(3),co(2):co(4)) = face;

subplot(1,2,1); imshow(faceIm);
title(sprintf('Subject %s Best Location: [%.2f %.2f] deg',subjConv(s).name,subjCov(s).bestDegX,subjCov(s).bestDegY));

subplot(1,2,2);
surface(X,Y,reshape(subjCov(s).result,length(sim.steps),length(sim.steps))); set(gca,'YDir','reverse');axis image;
title('All Sampled Locations');

niceSave(pwd,['/figures/' info.subjs{s} '_' ROI '_' fits.cond]);
save(['results/singleSubj_' ROI '_' fits.cond '.mat'],'subjCov','sim');
close all;
toc
end