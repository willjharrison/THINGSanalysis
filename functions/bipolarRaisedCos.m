function filter = bipolarRaisedCos(radDist,ctr,cut,logOrLinear)
% filter = raisedCos(radDist,ctr,cut,logOrLinear)
% 
% Creates a raised log (or linear) cosine filter, similar to that described
% by Peli (1990), and used by Chung et al. radDist is a matrix of radial
% distances relative to the centre of an image (or 2d fourier transform).
% ctr is the peak frequency of the filter.

if strcmp(logOrLinear,'log')
    filter = (1 + cos(pi*((log(radDist)-log(ctr))/(log(cut)-log(ctr)))))/2;
    filter(radDist < cut | radDist > ctr^2/cut) = 0;
else
    
    thisCut = cut - ctr;

    thisMat = rad2deg(circ_dist(deg2rad(radDist),deg2rad(ctr)));
    firstHalf = (1 + cos(pi*((thisMat)/(thisCut))))/2;
    firstHalf(abs(thisMat) > abs(thisCut)) = 0;
 
    thisMat = rad2deg(circ_dist(deg2rad(radDist),deg2rad(ctr + 180)));
    secondHalf = (1 + cos(pi*((thisMat)/(thisCut))))/2;
    secondHalf(abs(thisMat) > abs(thisCut)) = 0;
    
    filter = firstHalf + secondHalf;
    
end

