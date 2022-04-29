function C = lab2rgb_array(L,a,b)

if nargin==1
    a = L(:,2); b = L(:,3); L = L(:,1);
end

% Convert CIE L*a*b* to CIE XYZ
V(:,:,2) = ( L + 16 ) / 116;  % (Y/Yn)^(1/3)
V(:,:,1) = a / 500 + V(:,:,2);  % (X/Xn)^(1/3)
V(:,:,3) = V(:,:,2) - b / 200;  % (Z/Zn)^(1/3)

Z = V.^3; 

% Correction for small XYZ
Z(Z <= 0.008856) = (V(Z <= 0.008856) - 16 / 116) / 7.787;

% Adjust for white point (D65, CIE 2 Deg Standard Observer)
Zn = cat(3,95.047, 100.00, 108.883);
Z = Z.*repmat(Zn,size(Z,1),size(Z,2),1);

% Convert CIE XYZ to Rec 709 RGB
M = [ 3.2406  -1.5372  -0.4986;
     -0.9689   1.8758   0.0415;
      0.0557  -0.2040   1.0570];

Zadjust = reshape(Z,size(Z,1)*size(Z,2),3);

R = Zadjust*M' / 100;

% Correct for non-linear output of display 
Q = R*12.92;
Q(R > 0.0031308) = 1.055 * R(R > 0.0031308).^(1/2.4) - 0.055;

% Scale to range 1-255
C = round(Q*255); C(C>255)=255; C(C<0)=0;
C = reshape(C,size(Z,1),size(Z,2),3);
