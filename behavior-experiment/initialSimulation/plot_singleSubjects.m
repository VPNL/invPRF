ROI = 'pFus_faces';
conds = {'Inverted' 'Upright'};

for c = 1:length(conds)
load(['results/singleSubj_' ROI '_' conds{c} '.mat']);
load(pRFfile([raid 'invPRF/'],'fixPRF',20,'photo','kayCSS',{'rh' 'lh'},''));
res = roi(1).fits(1).res; ppd = roi(1).fits(1).ppd;

res = res+ppd*5; % give a little buffer around our image/make indexing easier

load('stims/internalAvg.mat'); if strcmp(conds{c},'Inverted'); avgFace = flipud(avgFace); end
faceSize = sim.faceSize*ppd;
face = imresize(avgFace,[faceSize faceSize]);

for s = 1:length(subjCov)
    faceIm = zeros(res+1,res+1); co = CenterRectOnPoint([1 1 faceSize faceSize],round(subjCov(s).bestDegY *ppd)+res/2+1,round(subjCov(s).bestDegX *ppd)+res/2+1);
    faceIm(co(1):co(3),co(2):co(4)) = face;
    bestIm(s,:,:) = imcomplement(faceIm);
end

subplot(2,2,1+2*(c-1));
imagesc(squeeze(mean(bestIm))); colormap('gray'); set(gca,'visible','off');
axis([0 res 0 res]);set(gca,'YDir','reverse'); hline(res/2+1,'b:'); vline(res/2+1,'b:'); hold on; axis image;

subplot(2,2,2+2*(c-1));
axis([-res/2/ppd res/2/ppd -res/2/ppd res/2/ppd]);set(gca,'YDir','reverse'); hline(0,'b:'); vline(0,'b:'); hold on; axis square;
s = scatterCent([subjCov.bestDegX],[subjCov.bestDegY],condColors(1),[],[],[],12,1,0);
text([subjCov.bestDegX]+.1, [subjCov.bestDegY]+.1, info.subjs', 'Fontsize', 5);


grid on; 
try
    xticks(-res/2/ppd:res/2/ppd); yticks(-res/2/ppd:res/2/ppd);
catch 
    set(gca,'XTick',-res/2/ppd:res/2/ppd);set(gca,'YTick',-res/2/ppd:res/2/ppd);
end
title(conds{c});

end
superTitle(['"Best" Face Position from ' ROI ' pRF Coverage, N = ' num2str(length(subjCov))],12,.95);
