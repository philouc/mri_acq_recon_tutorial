% function f = myFFT (x);
%
% computes orthogonal Fourier Transform
%
% developper : Nicolas Chauffert

function f= myFFT (x)
f=fftn(x) / sqrt(numel(x));