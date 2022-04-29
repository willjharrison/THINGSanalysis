function C = colouriseOris_withAmp(im,amp,Lmax)
% takes an input "im" that is orientations between -pi and pi and converts
% it to LAB coloursspace. This file will additiopnally take the amplitude
% of the colour, change the saturation. C is an imSize x imSize x 3 rgb
% image.

if nargin < 3
    Lmax = 70;
end

% convert to colour
L = amp./max(amp(:))*Lmax;
[a,b] = pol2cart(im,L);
C = lab2rgb_array(L,a,b)/255;