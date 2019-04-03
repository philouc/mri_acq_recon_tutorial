% function f = myIFFT (x);
%
% computes orthogonal Inverse Fourier Transform
%
% developper : Nicolas Chauffert

function f= myIFFT (x)
f=ifftn(x) * sqrt(numel(x));