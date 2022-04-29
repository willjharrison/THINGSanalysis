function [angDist,radDist] = polarDistFun(imSize)
% [angDist,radDist] = polarDistFun(imSize)
% 
% Calculates the radial and angular distances from the centre of a square
% image with size specified by imSize. angDist is in degrees, not radians.
% 
% Updated to allow different x/y Feb 2018

if length(imSize) == 1
    imSize(2) = imSize(1);
end

smallestSize = min(imSize);

[X,Y] = meshgrid(linspace(-smallestSize/2,smallestSize/2-1,imSize(2)),linspace(-smallestSize/2,smallestSize/2-1,imSize(1)));

% [X,Y]=meshgrid(-imSize(2)/2:imSize(2)/2-1,-imSize(1)/2:imSize(1)/2-1); % 2D matrix of radial distances from centre

radDist=(X.^2+Y.^2).^0.5; % radial distance from centre
angDist=rad2deg(atan2(-Y, X)); % orientation around image

