%% analysis file to test model but with normalised images
% this file will load every THINGS image, find it's RMS contrast, mean
% luminance, and orientation/SF bandlimited energy. These data are
% summarised, and also used to classify each image based on only these
% metrics. Different from the main analysis, however, images will be
% normalised in their luminance and contrast as described in the
% manuscript. Therefore, decodabiity based on these metrics should be
% vastly lower than from the main analysis. 
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
imStats.image_concept_index = image_concept_index;
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % normalise luminance and contrast
        im = im - mean(im(:));
        im = im/std(im(:))*.2; % 0.2 is arbitrary
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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

%% save variables as csv file
imStats.oriEn = imStats.oriEn.';
imStats.sfEn = imStats.sfEn.';
THINGS_imstats_normed = struct2table(imStats);
writetable(THINGS_imstats_normed, './output/THINGS_imstats_normed.csv');

%% COMPUTE MEANS INCLUDING LEAVE ONE OUTS

allLumMeans = qStats(imStats.lum,imStats.image_concept_index,6);
allRMSMeans = qStats(imStats.RMS,imStats.image_concept_index,6);

% n - 1 mean function (leave one out mu, looMu)
looMu = @(mu,x,n) (mu - x*(1/n))*(n./(n-1));

% n - 1 std function (leave one out SD, looSD)
looSD = @(sd,mu,x,n) sqrt(...
    (((1 - n).^2)*sd^2 - n*((mu - x).^2))...
    ./((2 - n).*(1 - n)));

looMu_lum = [];
looSD_lum = [];
looMu_RMS =[];
looSD_RMS =[];
for i = 1:1854
    
    % luminance
    x = imStats.lum(imStats.image_concept_index==i);
    mu = allLumMeans(i,2);
    sd = allLumMeans(i,5);
    n = allLumMeans(i,4);
    looMu_lum = [looMu_lum; looMu(mu,x,n)];
    looSD_lum = [looSD_lum; looSD(sd,mu,x,n)];    
    
    % contrast
    x = imStats.RMS(imStats.image_concept_index==i);
    mu = allRMSMeans(i,2);
    sd = allRMSMeans(i,5);
    n = allRMSMeans(i,4);
    looMu_RMS = [looMu_RMS; looMu(mu,x,n)];
    looSD_RMS = [looSD_RMS; looSD(sd,mu,x,n)];
    
end

%% FIND DISTANCES

% distances from LOO means
lumDiff_loo = (imStats.lum - looMu_lum)./looSD_lum;
rmsDiff_loo = (imStats.RMS - looMu_RMS)./looSD_RMS;
lumRMSDiff_loo = sqrt(lumDiff_loo.^2 + rmsDiff_loo.^2);

% distances from every other mean (we include the self-group mean in the
% calculation here, but we'll get rid of it later)
lumDiff = (imStats.lum - repmat(allLumMeans(:,2)',length(imStats.lum),1))...
    ./repmat(allLumMeans(:,5)',length(imStats.lum),1);
rmsDiff = (imStats.RMS - repmat(allRMSMeans(:,2)',length(imStats.lum),1))...
    ./repmat(allRMSMeans(:,5)',length(imStats.lum),1);
lumRMSDiff = sqrt(lumDiff.^2 + rmsDiff.^2);

% remove the within-group difference scores because we use the LOO differences
idx = sub2ind(size(lumRMSDiff),(1:length(lumRMSDiff))', imStats.image_concept_index);
lumRMSDiffReduced = lumRMSDiff;
lumRMSDiffReduced(idx) = NaN;

%%
% compare LOO mean with every other mean
pairwiseClassCorrect = lumRMSDiff_loo < lumRMSDiff;

% average across all concepts
[s,n] = grpstats(pairwiseClassCorrect, imStats.image_concept_index,{'sum','numel'});
propCorrect = s./(n);

% average across diaganols - this is needed because the different numbers
% of images in different categories results in propCorrect being
% non-symmetric
propCorrect = (propCorrect + propCorrect')/2;

mean(propCorrect(:))



