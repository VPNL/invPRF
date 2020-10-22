% simulates overlap of pRF coverage and mean face features at different
% locations

clear all; close all;


expt = 'fixPRF'; % assumes all subject in pRFset

saveFig = 1;

minR2 = 20;          % cutoff for vox selection
ROI= 'IOG_faces';

whichStim = 'photo';%'internal';%
whichModel = 'kayCSS';%'cssShift';%
hems = {'rh' 'lh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in data, sort plotVox           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;

% now we load in the data from both hemispheres, and threshold across
load(pRFfile('',expt,minR2,whichStim,whichModel,hems,''));
fits = roi(cellNum(ROI,info.ROIs)).fits; % one ROI at a time, for now
res = roi(1).fits(1).res; ppd = roi(1).fits(1).ppd;

res = res+ppd*5; % give a little buffer around our image/make indexing easier

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make coverage map                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y]=meshgrid(-res/2:res/2);

c=1; % inverted
for v = 1:length(fits(c).vox)
    cond(c).cov(v,:,:) = PRF(X,Y,fits(c).vox(v).XYdeg(1)*ppd,fits(c).vox(v).XYdeg(2)*ppd,fits(c).vox(v).size*ppd);
end

% to deal with the issue of plot/imagesc using different coordinate systems
cond(c).covIm = flipud(squeeze(mean(cond(c).cov))); 
subplot(1,3,1); imshow(cond(c).covIm); colormap(mrvColorMaps('hot'));title([ROI ' coverage']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation parameters                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sim.steps = -3:.1:3; [X,Y]=(meshgrid(sim.steps));
sim.degPos = [Y(:) X(:)];
sim.centers = sim.degPos*ppd+res/2+1; 
sim.faceSize = 3.2;

% load & resize face
load('internalAvg.mat'); if c == 1 avgFace = flipud(avgFace); end
faceSize = sim.faceSize*ppd;
face = imresize(avgFace,[faceSize faceSize]);

for n = 1:length(sim.centers)
% make current face im
faceIm = zeros(res+1,res+1); co = CenterRectOnPoint([1 1 faceSize faceSize],sim.centers(n,1),sim.centers(n,2));
faceIm(co(1):co(3),co(2):co(4)) = face;

% subplot(1,3,2); imshow(faceIm);

featCov = faceIm.*cond(c).covIm; %subplot(1,3,3); imshow(featCov);
sim.result(n) = sum(featCov(:));
end

best = find(sim.result==max(sim.result)); worst = find(sim.result==min(sim.result));

% make best face im
faceIm = zeros(res+1,res+1); co = CenterRectOnPoint([1 1 faceSize faceSize],sim.centers(best,1),sim.centers(best,2));
faceIm(co(1):co(3),co(2):co(4)) = face;

subplot(1,3,2); imshow(faceIm); 
title(sprintf('Best Location: [%.2f %.2f] deg',sim.degPos(best,2),sim.degPos(best,1)));

subplot(1,3,3);
surface(X,Y,reshape(sim.result,length(sim.steps),length(sim.steps))); set(gca,'YDir','reverse');axis image;
title('All Sampled Locations');

superTitle(['pRF coverage x face features (' num2str(sim.faceSize) 'deg)'],12,.05); if onLaptop playSound; end