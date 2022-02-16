clear
close all
clc

load('THINGS_imstats.mat');
nCats = length(unique(imStats.image_concept_index));

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

% pull out only one half of symmetric matrix
idx = logical(1 - triu(ones(size(propCorrect))));
halfProp = propCorrect(idx);

%% organise for plotting
p = propCorrect(:);
concept = numberList(nCats,nCats);
propForDgraph = table(concept, p);
writetable(propForDgraph,'propForDgraph_revised_2.csv') 

% corrdinates for plotting
[x,y] = ind2sub(nCats, 1:nCats^2);
x = x';
y = y';
coords = table(x, y);
writetable(coords, 'coords_revised_2.csv');

% save matrix as an image with colour bar
bigMatrix = reshape(p,nCats,nCats);
lumMatrix = abs(bigMatrix*2 - 1).^(1/3);
colMatrix = round(bigMatrix)*pi - pi/2;
imwrite(flipud(colouriseOris_withAmp(colMatrix,lumMatrix,100)),'colIm_revised_2.png');

%% col bar
colBar = round(0:.001*2:1);
lumBar = (abs(linspace(0,1,length(colBar))*2 - 1)).^(1/3);

colBarIm = colouriseOris_withAmp(kron(colBar*pi - pi/2,ones(50,1)), kron(lumBar,ones(50,1)),100);
imwrite(colBarIm, 'colBarIm_revised_2.png');

%% zoom in on section of matrix
% smaller section
w = 101:120;
smallP = propCorrect(w,w);

smallP = smallP(:);
smallConcept = min(w) - 1 + numberList(length(w), length(w));
smallPropForDgraph = table(smallConcept, smallP);
writetable(smallPropForDgraph,'smallPropForDgraph_revised_2.csv') 

smallMatrix = reshape(smallP,20,20);
imwrite(flipud(colouriseOris_withAmp(...
    kron(round(smallMatrix)*pi - pi/2,ones(2)),...
    kron(abs(smallMatrix*2 - 1).^(1/3),ones(2)),...
    100)),'smallMatrixIm_revised_2.png');

load('/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Variables/unique_id.mat')
labels = string(unique_id(w));
catNames = string(unique_id);

%% which images are easiest/hardest to classify?
load('/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Variables/image_paths.mat')
thingsDir = '/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Main/';

nBestWorst = 5; % arbitrary

% most distinguished images
meanImAcc = mean(pairwiseClassCorrect,2);
[maxAcc,mostAccIdx] = maxk(meanImAcc,nBestWorst);
mostAcc = image_paths(mostAccIdx);

% least distinguished images
[minAcc,leastAccIdx] = mink(meanImAcc,nBestWorst);
leastAcc = image_paths(leastAccIdx);

% load the images and arrange them sensibly
newImSize = 128;
ims = zeros(newImSize,newImSize,nBestWorst,2);
 for imLoop = 1:nBestWorst
        
        whichIm = imLoop;
        
        ims(:,:,imLoop,1) = imresize((rgb2gray(double(imread([thingsDir mostAcc{whichIm}]))/255).^(1))*2 - 1,[newImSize newImSize]);
        ims(:,:,imLoop,2) = imresize((rgb2gray(double(imread([thingsDir leastAcc{whichIm}]))/255).^(1))*2 - 1,[newImSize newImSize]);
                
 end
ims(ims(:)>1) = 1;
ims(ims(:)<-1) = -1;
imwrite(unpackIm(sqrt(ims/2 + .5),0,0),'bestWorstIms_revised.png')