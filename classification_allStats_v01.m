clear
close all
clc

load('./output/THINGS_imstats.mat');
nCats = length(unique(imStats.image_concept_index));

%% COMPUTE MEANS INCLUDING LEAVE ONE OUTS

allLumMeans = qStats(imStats.lum,imStats.image_concept_index,6);
allRMSMeans = qStats(imStats.RMS,imStats.image_concept_index,6);
allSFStats = qStats([imStats.sfEn(:,1:8)], imStats.image_concept_index,6);
allOriStats = qStats([imStats.oriEn], imStats.image_concept_index,6);

% n - 1 mean function (leave one out mu, looMu)
looMu = @(mu,x,n) (mu - x.*(1./n)).*(n./(n-1));

% n - 1 std function (leave one out SD, looSD)
looSD = @(sd,mu,x,n) sqrt(...
    (((1 - n).^2).*sd.^2 - n.*((mu - x).^2))...
    ./((2 - n).*(1 - n)));

% stfu :)
looMu_lum = [];
looSD_lum = [];
looMu_RMS =[];
looSD_RMS =[];
looMu_SFs = [];
looSD_SFs = [];
looMu_Ori = [];
looSD_Ori = [];
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

    % all SFs
    x = imStats.sfEn(imStats.image_concept_index==i,1:8); % only 8 Sfs
    mu = allSFStats(i,2:9);
    sd = allSFStats(i,26:33);
    n = allSFStats(i,18:25);
    looMu_SFs = [looMu_SFs; looMu(mu,x,n)];
    looSD_SFs = [looSD_SFs; looSD(sd,mu,x,n)];

    % all RMS
    x = imStats.oriEn(imStats.image_concept_index==i,:);
    mu = allOriStats(i,2:18);
    sd = allOriStats(i,53:69);
    n = allOriStats(i,18:34);
    looMu_Ori = [looMu_Ori; looMu(mu,x,n)];
    looSD_Ori = [looSD_Ori; looSD(sd,mu,x,n)];

end

%% FIND DISTANCES

% distances from LOO means
lumDiff_loo = (imStats.lum - looMu_lum)./looSD_lum;
rmsDiff_loo = (imStats.RMS - looMu_RMS)./looSD_RMS;
sfDiff_loo = (imStats.sfEn(:,1:8) - looMu_SFs)./looSD_SFs;
oriDiff_loo = (imStats.oriEn - looMu_Ori)./looSD_Ori;

allDiff_loo = sqrt(lumDiff_loo.^2 + rmsDiff_loo.^2 + ...
    sum(sfDiff_loo.^2,2) + sum(oriDiff_loo.^2,2));

%%
% distances from every other mean (we include the self-group mean in the
% calculation here, but we'll get rid of it later)
lumDiff = (imStats.lum - repmat(allLumMeans(:,2)',length(imStats.lum),1))...
    ./repmat(allLumMeans(:,5)',length(imStats.lum),1);
rmsDiff = (imStats.RMS - repmat(allRMSMeans(:,2)',length(imStats.lum),1))...
    ./repmat(allRMSMeans(:,5)',length(imStats.lum),1);

% same as above, but in larger array which requires a loop
SFs = 2.^(0:8); % labels for use later
sfDiff = zeros(length(lumDiff), size(lumDiff,2), length(SFs)-1);
for i = 1:length(SFs)-1

    % scale by linear distance
%     sfDiff(:,:,i) = (imStats.sfEn(:,i) - repmat(sfOriStats(:,1+i)',length(imStats.sfEn),1));

    % alternatively, distances can be scaled by the SD of group data
    sfDiff(:,:,i) = (imStats.sfEn(:,i) - repmat(allSFStats(:,1+i)',length(imStats.sfEn),1))...
        ./repmat(allSFStats(:,25+i)',length(imStats.sfEn),1);
     
end

% same as above, but in larger array which requires a loop
oris = 0:22.5/2:(180);
oriDiff = zeros(length(lumDiff), size(lumDiff,2), length(oris)-1);
for i = 1:length(oris)-1

    % scale by linear distance
%     oriDiff(:,:,i) = (imStats.oriEn(:,i) - repmat(sfOriStats(:,10 + i)',length(imStats.sfEn),1));
    
    % alternatively, distances can be scaled by the SD of group data
    oriDiff(:,:,i) = (imStats.oriEn(:,i) - repmat(allOriStats(:,1 + i)',length(imStats.sfEn),1))...
        ./repmat(allOriStats(:,52 + i)',length(imStats.sfEn),1);
    
end

%%

allDiff = sqrt(lumDiff.^2 + rmsDiff.^2 + sum(sfDiff.^2,3) + sum(oriDiff.^2,3));

% remove the within-group difference scores and use the LOO differences
idx = sub2ind(size(allDiff),(1:length(allDiff))', imStats.image_concept_index);
allDiffReduced = allDiff;
allDiffReduced(idx) = allDiff_loo; % note this line is different from the pairwise script

[~,minDiffIdx] = min(allDiffReduced,[],2);
correctClassifications_allStats = minDiffIdx == imStats.image_concept_index;

[sum(correctClassifications_allStats) mean(correctClassifications_allStats)]

chance = 1/1854;
pCor_allStats = 1 - binocdf(sum(correctClassifications_allStats),length(correctClassifications_allStats),chance);
