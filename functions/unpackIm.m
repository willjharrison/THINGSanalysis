function imOut = unpackIm(im,showIm, normContrast)
% imOut = unpackIm(im,showIm, normContrast)
% 
% For a multi-dimensional array of images (up to 4D), this will turn the
% image into 2D laid out sensibly. normContrast set to 1 normalises the
% contrast range for each panel, and puts the entire image on the scale 0 -
% 1.

if nargin == 1
    showIm = 0;
end

if nargin < 3
    normContrast = 0;
end

imSize = size(im);

dim = length(imSize);

% normalise contrast
if normContrast
    if dim == 3
        for i = 1:imSize(dim)
            im(:,:,i) = maxContrast(im(:,:,i))/2 + .5;
        end
    elseif dim == 4
        for i = 1:imSize(3)
            for j = 1:imSize(4)
                im(:,:,i,j) = maxContrast(im(:,:,i,j))/2 + .5;
            end
        end
    end
end

% reshape image
if dim == 3
    imOut = reshape(im,imSize(1),imSize(2)*imSize(3));
elseif dim == 4
    im = permute(im,[1 3 2 4]);
    imOut = reshape(im,imSize(2)*imSize(3),imSize(1)*imSize(4));
else
    fprintf('\rYou tit! Image dimensions not appropriate for this function!\r')
    imOut=im;
end

if showIm == 1 && normContrast == 0
    imshow(imOut,[]);
elseif showIm == 1 && normContrast ~= 0
    imshow(imOut);
end