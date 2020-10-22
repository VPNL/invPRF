clear all; close all;

which = {'full' 'outline' 'internal' 'eyes'};


%%% for some basic params
load internalAvg.mat
faceSize = size(avgFace,1); ppd = faceSize/3.2;
imCenter = faceSize/2;

for n = 1:length(which)
    im = imread ([dirOf(pwd) 'stimuli/' which{n} '/face7.png']);
    thresh = imread(['imOffset_mockup/' which{n} '_thresh.png']);
    thresh = thresh./255;
    subplot(2,4,n);

imshow(im); hold on;
[centX,centY] = weightedCent(im,thresh);
h = hline(centY,'b-');set(h,'LineWidth',2)

offset(n) = (centY-imCenter)/ppd;
title(sprintf('%s: %.3fdeg from horiz.',which{n},offset(n)));

subplot(2,4,n+4); imshow(flipud(im));hold on;
[centX,centY] = weightedCent(flipud(im),flipud(thresh)); 
h = hline(centY,'b-');set(h,'LineWidth',2)
end
