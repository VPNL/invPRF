% simulates overlap of pRF coverage and mean face features at different
% locations

%clear all; close all;

expt = 'fixPRF'; % assumes all subject in pRFset
minR2 = 20;          % cutoff for vox selection
ROI= 'mFus_faces';

whichStim = 'photo';%'outline';%'internal';%
whichModel = 'kayCSS';%'cssShift';%
hems = {'rh' 'lh'};
c=1 ; % 1 = inverted

load(['results/' ROI '_c' num2str(c) '_sim.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in data, sort plotVox           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;

% now we load in the data from both hemispheres, and threshold across
load(pRFfile([raid 'invPRF/'],expt,minR2,whichStim,whichModel,hems,''));
fits = roi(cellNum(ROI,info.ROIs)).fits; % one ROI at a time, for now
res = roi(1).fits(1).res; sim.ppd = roi(1).fits(1).ppd;

sim.res = res+sim.ppd*1; % give a little buffer around our image/make indexing easier

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation parameters                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sim.steps = -3:.1:3; [X,Y]=(meshgrid(sim.steps));
sim.degPos = [Y(:) X(:)];
sim.centers = sim.degPos*sim.ppd+sim.res/2+1;
sim.faceSize = 3.2;

% load & resize face
load('stims/internalAvg.mat'); if c==1 avgFace = flipud(avgFace); end
faceSize = sim.faceSize*sim.ppd;
face = imresize(avgFace,[faceSize faceSize]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make coverage map                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% [XX,YY]=meshgrid(-sim.res/2:sim.res/2);
% 
% for s = 1%:sim.numSims
%     tic
%     draw(s).voxInd = randsample(length(fits(1).vox),sim.drawVox);
%     for vv = 1:sim.drawVox
%         v = draw(s).voxInd(vv);
%         cov(vv,:,:) = PRF(XX,YY,fits(c).vox(v).XYdeg(1)*sim.ppd,fits(c).vox(v).XYdeg(2)*sim.ppd,fits(c).vox(v).size*sim.ppd);
%     end
%     
%     % to deal with the issue of plot/imagesc using different coordinate systems
%     covIm = flipud(squeeze(mean(cov)));
%     subplot(1,3,1); imshow(covIm); colormap(mrvColorMaps('hot'));title([ROI ' coverage']);
%     
%     for n = 1:length(sim.centers)
%         % make current face im
%         faceIm = zeros(sim.res+1,sim.res+1); co = CenterRectOnPoint([1 1 faceSize faceSize],sim.centers(n,1),sim.centers(n,2));
%         faceIm(co(1):co(3),co(2):co(4)) = face;
%         
%          subplot(1,3,2); imshow(faceIm);
%         
%         featCov = faceIm.*covIm; subplot(1,3,3); imshow(featCov);
%     end
% end
% if onLaptop playSound; end
% subplot(1,3,3)
[bestPos,sim] = drawBestFace(draw,sim,face);
superTitle(ROI,12,.05)
niceSave(pwd,['/figures/' ROI '_' whichStim '_' fits(c).cond],[],[],{'svg'});
