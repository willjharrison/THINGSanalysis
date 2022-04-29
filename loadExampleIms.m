%% load some example THINGS images
% images are shown in colour, greyscale, and as an FFT amplitude

clear
close all
clc

rSeed = rng(1); % for reproducibility

%% setup filters
stimuli.imSize = 256;

%% find all images
encodingGamma = 1;
load('./THINGS/Variables/image_concept_index.mat')
load('./THINGS/Variables/image_paths.mat')
load('./THINGS/Variables/unique_id.mat')
thingsDir = './THINGS/Main/';

nCats = length(unique(image_concept_index));

% randomly select some categories
randCats = randperm(nCats,10);

%% grab ims
im = zeros(stimuli.imSize,stimuli.imSize,3,length(randCats));
greyIm = zeros(stimuli.imSize,stimuli.imSize,length(randCats));
imFFT = greyIm;
counter = 1;

for catLoop = 1:length(randCats)
    
    whichCat = randCats(catLoop);
    allTargIms = image_paths(image_concept_index == whichCat);
    nIms = length(allTargIms);
    
    % just grab first im
    whichIm = 1;
    
    im(:,:,:,catLoop) = imresize((double(imread([thingsDir allTargIms{whichIm}])))/255,[stimuli.imSize stimuli.imSize]);
    
    greyIm(:,:,catLoop) = rgb2gray(im(:,:,:,catLoop));
    
    imFFT(:,:,catLoop) = log(abs(fftshift(fft2(greyIm(:,:,catLoop)))));
end

% fix possible errors from resize and other calcs
im(im(:)>1) = 1;
im(im(:)<0) = 0;


originalIms = reshape(permute(im,[1 2 4 3]), stimuli.imSize, stimuli.imSize*length(randCats),3);
normFFT = imFFT;
normFFT = normFFT - min(normFFT(:));
normFFT = normFFT/max(normFFT(:));
allIms = [originalIms; repmat(unpackIm(greyIm,0,0),1,1,3); repmat(unpackIm(normFFT,0,0),1,1,3)];
imshow(allIms)

imwrite(allIms, './output/allIms.png')