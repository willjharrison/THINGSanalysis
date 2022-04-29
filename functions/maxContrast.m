function newIm = maxContrast(im,newRange)
% im has to be in the range of -1 to 1

if nargin < 2
    newRange = 0;
end

    
    maxL = max(im(:));
    minL = min(im(:));
    
    if maxL > abs(minL)
        
        newIm = im/maxL;
        
    else
        
        newIm = im/abs(minL);
        
    end
if newRange == 1
    newIm = newIm - min(newIm(:));
    newIm = newIm/max(newIm(:));
end