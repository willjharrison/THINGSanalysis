%% main analysis file
% this file will load every THINGS image, find it's RMS contrast, mean
% luminance, and orientation/SF bandlimited energy. 
% 
% Please get in touch if you have any questions or find any errors:
% willjharri@gmail.com
% 
% If you'd like me to do various image analyses as part of your project,
% I'd be more than happy to collaborate!

clear
close all
clc

%% setup filters
% set up a single image size for all images
stimuli.imSize = 512;

% create circular aperture
[angDist, radDist] = polarDistFun(stimuli.imSize); % this function finds the angle and distance of every point
mask = ...
    createLowPassFilter(radDist,stimuli.imSize/2 - 6*2,6*2);

%% filters to calculate orientation energy

% fully wrap orientation for plotting
oris = 0:22.5/2:(180);
oriBw = 22.5/2;
for i = 1:length(oris)
    
    oriFilts(:,:,i) = fftshift(bipolarRaisedCos(angDist, oris(i), oris(i) - oriBw, 'linear'));
    
end

% select log-spaced SF filters
SFs = 2.^(0:8); % labels for use later

% transform distances to log steps for simplicity
sradDist = floor(log2(fftshift(radDist)));
sradDist(1) = 0; % ignore DC of images

%% find all THINGS images
encodingGamma = 2;
load('./THINGS/Variables/image_concept_index.mat')
load('./THINGS/Variables/image_paths.mat')
load('./THINGS/Variables/unique_id.mat')
thingsDir = './THINGS/Main/';

nCats = length(unique(image_concept_index));

%% grab memory for outputs
imStats.oriEn = zeros(length(oris), length(image_concept_index));
imStats.sfEn = zeros(length(SFs), length(image_concept_index));
imStats.RMS = zeros(length(image_concept_index),1);
imStats.lum = zeros(length(image_concept_index),1);

%% do calcs
counter = 1;
tic
for catLoop = 1:nCats
    
    whichCat = catLoop;
    allTargIms = image_paths(image_concept_index == whichCat);
    nIms = length(allTargIms);
    
    for imLoop = 1:nIms
        
        whichIm = imLoop;
        
        % load im, greyscale, gamma correct, and set in range -1 to 1
        im = imresize((rgb2gray(double(imread([thingsDir allTargIms{whichIm}]))/255).^(encodingGamma))*2 - 1,[stimuli.imSize stimuli.imSize]);
        
        % clip bad values from imresize
        im(im(:)>1) = 1;
        im(im(:)<-1) = -1;
        
        % RMS contrast
        imStats.RMS(counter,1) = std(im(:));
        
        % mean luminance, relative to mid-grey
        imStats.lum(counter,1) = mean(im(:));
        
        % mask im for other calcs
        maskedIm = im.*mask;
        
        % find orientation energy
        imFFT = abs(fft2(maskedIm)).*fftshift(mask); % repeating the mask gets rid of high SF obliques
        imFFT(1) = 0;
        thisOriEn = oriFilts.*imFFT;
        imStats.oriEn(:,counter) = squeeze(sum(sum(thisOriEn)));
        
        % find SF energy
        thisSFEn = grpstats(imFFT(:), sradDist(:));
        imStats.sfEn(:,counter) = thisSFEn;
                
        counter = counter + 1;
    end
    
    % update me please
    if mod(catLoop,10)
        
        [catLoop mod(toc,60)]
        
    end
end

% save variables as csv file
% imStats.oriEn = imStats.oriEn.' 
% imStats.sfEn = imStats.sfEn.'
THINGS_imstats = struct2table(imStats);
mkdir('./output')
writetable(THINGS_imstats, './output/THINGS_imstats.csv');

%% summarise main data
% This section was removed during revisions and replaced by other files
% named "classification..."
%  
% 
% meanStats = qStats([imStats.lum imStats.RMS], imStats.image_concept_index,4);
% sfOriStats = qStats([imStats.sfEn imStats.oriEn], imStats.image_concept_index,4);
% 
% % spit out in format useful for datagraph
% sfMeans = [repelem(SFs(1:end-1)',nCats) repmat((1:nCats)',8,1) reshape(sfOriStats(:,2:9),numel(sfOriStats(:,2:9)),1)];
% oriMeans = [repelem(oris',nCats)-90 repmat((1:nCats)',length(oris),1) reshape(sfOriStats(:,11:27),numel(sfOriStats(:,11:27)),1)];
% 
% %% classify each image according to nearest category mean luminance, contrast, SF, and orientation energy
% 
% tic;
% 
% % find difference of each value from the mean of each category
% % scale by linear distance
% % lumDiff = imStats.lum - repmat(meanStats(:,2)',length(imStats.lum),1);
% % rmsDiff = imStats.RMS - repmat(meanStats(:,3)',length(imStats.lum),1);
% 
% % alternatively, distances can be scaled by the SD of group data
% lumDiff = (imStats.lum - repmat(meanStats(:,2)',length(imStats.lum),1))...
%     ./repmat(meanStats(:,6)',length(imStats.lum),1);
% rmsDiff = (imStats.RMS - repmat(meanStats(:,3)',length(imStats.lum),1))...
%     ./repmat(meanStats(:,7)',length(imStats.lum),1);
% 
% % same as above, but in larger array which requires a loop
% sfDiff = zeros(length(lumDiff), size(lumDiff,2), length(SFs)-1);
% for i = 1:length(SFs)-1
% 
%     % scale by linear distance
% %     sfDiff(:,:,i) = (imStats.sfEn(:,i) - repmat(sfOriStats(:,1+i)',length(imStats.sfEn),1));
% 
%     % alternatively, distances can be scaled by the SD of group data
%     sfDiff(:,:,i) = (imStats.sfEn(:,i) - repmat(sfOriStats(:,1+i)',length(imStats.sfEn),1))...
%         ./repmat(sfOriStats(:,53+i)',length(imStats.sfEn),1);
%      
% end
% 
% % same as above, but in larger array which requires a loop
% oriDiff = zeros(length(lumDiff), size(lumDiff,2), length(oris)-1);
% for i = 1:length(oris)-1
% 
%     % scale by linear distance
% %     oriDiff(:,:,i) = (imStats.oriEn(:,i) - repmat(sfOriStats(:,10 + i)',length(imStats.sfEn),1));
%     
%     % alternatively, distances can be scaled by the SD of group data
%     oriDiff(:,:,i) = (imStats.oriEn(:,i) - repmat(sfOriStats(:,10 + i)',length(imStats.sfEn),1))...
%         ./repmat(sfOriStats(:,62 + i)',length(imStats.sfEn),1);
%     
% end
% 
% toc
% 
% % do classification (find linear distances of each image from concept
% % means)
% % first, try just lum and contrast distance:
% totalDiff_lumCont = sqrt(lumDiff.^2 + rmsDiff.^2);
% 
% [~,minDiffIdx] = min(totalDiff_lumCont,[],2);
% correctClassifications_lumCont = minDiffIdx == imStats.image_concept_index;
% 
% [sum(correctClassifications_lumCont) mean(correctClassifications_lumCont)]
% 
% chance = 1/1854;
% pCor_lumCont = 1 - binocdf(sum(correctClassifications_lumCont),length(correctClassifications_lumCont),chance);
% 
% % then use all data
% totalDiff = sqrt(lumDiff.^2 + rmsDiff.^2 + sum(sfDiff.^2,3) + sum(oriDiff.^2,3));
% 
% [~,minDiffIdx] = min(totalDiff,[],2);
% correctClassifications = minDiffIdx == imStats.image_concept_index;
% 
% [sum(correctClassifications) mean(correctClassifications)]
% pCor_lumContAll = 1 - binocdf(sum(correctClassifications),length(correctClassifications),chance);
% 
% %% try classifying which of two images using just luminance and contrast
% 
% % linear distances
% lumRMSDiff = sqrt(lumDiff.^2 + rmsDiff.^2);
% idx = sub2ind(size(lumRMSDiff),(1:length(lumRMSDiff))', imStats.image_concept_index);
% targDist = lumRMSDiff(idx);
% pairwiseClassCorrect = lumRMSDiff(idx) < lumRMSDiff;
% 
% % average across all concepts
% [s,n] = grpstats(pairwiseClassCorrect, imStats.image_concept_index,{'sum','numel'});
% propCorrect = s./(n);
% 
% % average across diaganols - this is needed because the different numbers
% % of images in different categories results in propCorrect being
% % non-symmetric
% propCorrect = (propCorrect + propCorrect')/2;
% 
% % organise for plotting
% p = propCorrect(:);
% concept = numberList(nCats,nCats);
% propForDgraph = table(concept, p);
% writetable(propForDgraph,'propForDgraph.csv') 
% 
% % corrdinates for plotting
% [x,y] = ind2sub(nCats, 1:nCats^2);
% x = x';
% y = y';
% coords = table(x, y);
% writetable(coords, 'coords.csv');
% 
% % save matrix as an image with colour bar
% bigMatrix = reshape(p,nCats,nCats);
% imwrite(flipud(colouriseOris_withAmp(.5 + bigMatrix/2,bigMatrix,100)),'colIm.png');
% 
% % col bar
% colBar = .5:.001:1;
% lumBar = linspace(0,1,length(colBar));
% 
% colBarIm = colouriseOris_withAmp(kron(colBar,ones(50,1)), kron(lumBar,ones(50,1)),100);
% imwrite(colBarIm, 'colBarIm.png');
% 
% %% zoom in on section of matrix
% % smaller section
% w = 101:120;
% smallP = propCorrect(w,w);
% 
% smallP = smallP(:);
% smallConcept = min(w) - 1 + numberList(length(w), length(w));
% smallPropForDgraph = table(smallConcept, smallP);
% writetable(smallPropForDgraph,'smallPropForDgraph.csv') 
% 
% smallMatrix = reshape(smallP,20,20);
% imwrite(flipud(colouriseOris_withAmp(.5 + kron(smallMatrix,ones(2))/2,kron(smallMatrix,ones(2)),100)),'smallMatrixIm.png');
% 
% labels = string(unique_id(w));
% catNames = string(unique_id);
% 
% %% which images are easiest/hardest to classify?
% 
% nBestWorst = 5; % arbitrary
% 
% % most distinguished images
% meanImAcc = mean(pairwiseClassCorrect,2);
% [maxAcc,mostAccIdx] = maxk(meanImAcc,nBestWorst);
% mostAcc = image_paths(mostAccIdx);
% 
% % least distinguished images
% [minAcc,leastAccIdx] = mink(meanImAcc,nBestWorst);
% leastAcc = image_paths(leastAccIdx);
% 
% % load the images and arrange them sensibly
% newImSize = 128;
% ims = zeros(newImSize,newImSize,nBestWorst,2);
%  for imLoop = 1:nBestWorst
%         
%         whichIm = imLoop;
%         
%         ims(:,:,imLoop,1) = imresize((rgb2gray(double(imread([thingsDir mostAcc{whichIm}]))/255).^(1))*2 - 1,[newImSize newImSize]);
%         ims(:,:,imLoop,2) = imresize((rgb2gray(double(imread([thingsDir leastAcc{whichIm}]))/255).^(1))*2 - 1,[newImSize newImSize]);
%                 
%  end
% ims(ims(:)>1) = 1;
% ims(ims(:)<-1) = -1;
% imwrite(unpackIm(sqrt(ims/2 + .5),0,0),'bestWorstIms.png')

%% calculate the 1/f slope of all images

% "simple" 1/f model
sfEn = imStats.sfEn(:,1:8); % drop 9th column because it's empty
oneOverFModel = @(cl,sf,alpha) cl./(sf.^alpha);
costFun = @(x) sum((sfEn(:) - oneOverFModel(x(1), repelem(SFs(1:8)', 26107,1),x(2))).^2);
startVals = [mean(imStats.sfEn(:,1)), 1];
[oneOverModelParams,imfval_sf] = fminsearch(costFun,startVals);

% GLMM:
sfEn = imStats.sfEn(:,1:8);
y = log(sfEn(:));
x = repelem(log(SFs(1:8)), length(imStats.sfEn))';
g = repmat(image_concept_index,8,1);

dat = table(x,y,g);

% note that the following linear formula can be arranged as y = b0/(x^(-b1)), after re-adjusting log transforms above
sfModel = fitglme(dat,'y ~ x + (x | g)'); 
sfFixef = sfModel.fixedEffects
sfRanef = sfModel.randomEffects;
rangeOfExps = sfFixef(2) + sfRanef(2:2:end);

%% model the oriented contrast of all images

oriEn = imStats.oriEn(:,1:end-1);
mOris = oris(1:end-1) - 90; % drop circular orientation and get out of Fourier space
oriModel = @(b,beta,cv,rho,theta) ...
    (cv - rho*abs(sin(theta*pi/180*2)).^b) .* (((beta-1)/2+1) + (beta-1)/2*(cos(theta*pi/180*2)));
costFun = @(x) sum((oriEn(:) - oriModel(x(1),x(2),x(3),x(4),repelem(mOris', 26107,1))).^2);
startVals = [1 1 mean(imStats.oriEn(:,1)), 10e5];
[oriModelParams,imfval_ori] = fminsearch(costFun,startVals);
oriModelParams(1) % b
oriModelParams(2) % beta
oriModelParams(3) % cv
oriModelParams(4) % rho

%% calculate RMS vs luminance correlations

y = imStats.RMS;
x = imStats.lum;
g = image_concept_index;

dat = table(x,y,g);

mdl = fitglme(dat,'y ~ x + x^2 + (x + x^2 | g)');

ranef = mdl.randomEffects;




