function lpFilter = createLowPassFilter(radDist,edge,halfPeriod) 
% Low pass filter with a cosine edge.
% Creates a circular area with radius of size halfPeriod, and cosine edges.

lpFilter = (1 + cos(pi*((radDist-edge)/(halfPeriod))))/2;

lpFilter(radDist<edge) = 1;
lpFilter(radDist>edge+halfPeriod) = 0;

% imshow(filter,[])