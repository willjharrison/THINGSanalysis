%% normalising luminance and RMS contrast example
% this code gives an example of how to normalise luminance and contrast
% across images. The basic method is to convert the image to LUV format
% (rather than RGB), subtract the luminance mean, and then divide by the
% RMS contrast, and then transform the LUV back to RGB.

clear
close all
clc

rSeed = rng(32); % for reproducibility

%% setup filters
stimuli.imSize = 256;

%% find all images
encodingGamma = 1; % don't use gamma correction because that is taken into consideration in the LUV conversion
load('./THINGS/Variables/image_concept_index.mat')
load('./THINGS/Variables/image_paths.mat')
load('./THINGS/Variables/unique_id.mat')
thingsDir = './THINGS/Main/';

nCats = length(unique(image_concept_index));

% randomly select some categories
randCats = randperm(nCats,5);

%% grab ims
im = zeros(stimuli.imSize,stimuli.imSize,3,length(randCats));
counter = 1;
tic
for catLoop = 1:length(randCats)
    
    whichCat = randCats(catLoop);
    allTargIms = image_paths(image_concept_index == whichCat);
    nIms = length(allTargIms);
    
    % just grab first im
    whichIm = 1;
    
    im(:,:,:,catLoop) = imresize((double(imread([thingsDir allTargIms{whichIm}])))/255,[stimuli.imSize stimuli.imSize]);
    
end

% fix possible errors from resize
im(im(:)>1) = 1;
im(im(:)<0) = 0;

%% normalise contrast
rms = .23; % arbitrary
normIms = zeros(stimuli.imSize,stimuli.imSize,3,length(randCats));
for imLoop = 1:length(randCats)
    
    % convert to YUV
    thisIm = (im(:,:,:,imLoop));
    
    Y = 0.299 * thisIm(:,:,1) + 0.587 * thisIm(:,:,2) + 0.114 * thisIm(:,:,3);
    U = -0.14713 * thisIm(:,:,1) - 0.28886 * thisIm(:,:,2) + 0.436 * thisIm(:,:,3);
    V = 0.615 * thisIm(:,:,1) - 0.51499 * thisIm(:,:,2) - 0.10001 * thisIm(:,:,3);
    thisIm = cat(3,Y,U,V);
    preIm = thisIm(:,:,1);
    
    % normalise contrast
    lumIm = thisIm(:,:,1); % grab luminance plane
    lumIm = (lumIm.^encodingGamma)*2 - 1; % gamma correct and put in range of -1 to 1
    
    % grab luminance and contrast while we're here
    % RMS contrast
    preRMS(imLoop,1) = std(lumIm(:));
        
    % mean luminance, relative to mid-grey
    preLum(imLoop,1) = mean(lumIm(:));

    lumIm = lumIm - mean(lumIm(:)); % remove mean
    lumIm = lumIm/std(lumIm(:)); % standardise RMS
    lumIm = lumIm*rms; % set RMS
    lumIm(lumIm(:)>1) = 1;
    lumIm(lumIm(:)<-1) = -1;
    
    lumIm = (lumIm/2 + .5).^(1/encodingGamma);
    
    % return to image
    thisIm(:,:,1) = lumIm;
        
    % convert to rgb
    normIms(:,:,:,imLoop) = yuv2rgb(thisIm);
    
    % check to see if normalisation worked
    checkIm = (rgb2gray(normIms(:,:,:,imLoop)).^encodingGamma)*2 - 1;
        
    % RMS contrast
    checkRMS(imLoop,1) = std(checkIm(:));
        
    % mean luminance, relative to mid-grey
    checkLum(imLoop,1) = mean(checkIm(:));

end

% just check that we did what we wanted.
imInfo = table(preLum,preRMS,checkLum,checkRMS);

% save all that stuff.
originalIms = reshape(permute(im,[1 2 4 3]), stimuli.imSize, stimuli.imSize*length(randCats),3);
finalIms = reshape(permute(normIms,[1 2 4 3]), stimuli.imSize, stimuli.imSize*length(randCats),3);
imwrite(originalIms,'./output/originalIms.png')
imwrite(finalIms,'./output/finalIms.png')
imwrite(im(:,:,:,imLoop),'./output/exampleImCol.png')
imwrite(preIm,'./output/exampleImY.png')
imwrite(lumIm,'./output/exampleNormY.png')
imwrite(maxContrast(thisIm(:,:,2))/2 + .5,'./output/exampleImU.png')
imwrite(maxContrast(thisIm(:,:,3))/2 + .5,'./output/exampleImV.png')
figure;
imshow([originalIms;finalIms])