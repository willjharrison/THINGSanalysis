clear
close all
clc

load('./output/THINGS_imstats.mat');
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

% remove the within-group difference scores and use the LOO differences
idx = sub2ind(size(lumRMSDiff),(1:length(lumRMSDiff))', imStats.image_concept_index);
lumRMSDiffReduced = lumRMSDiff;
lumRMSDiffReduced(idx) = lumRMSDiff_loo; % note this line is different from the pairwise script

[~,minDiffIdx] = min(lumRMSDiffReduced,[],2);
correctClassifications_lumCont = minDiffIdx == imStats.image_concept_index;

[sum(correctClassifications_lumCont) mean(correctClassifications_lumCont)]

chance = 1/1854;
pCor_lumCont = 1 - binocdf(sum(correctClassifications_lumCont),length(correctClassifications_lumCont),chance);

