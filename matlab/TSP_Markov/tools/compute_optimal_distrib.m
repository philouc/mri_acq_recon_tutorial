function [pi L]=compute_optimal_distrib(N,options,level)
%function [pi L]=compute_optimal_distrib(N,options,level)
%
% Computes the optimal distribution on a cuboid of size N = [n1 n2 ...]
% N=size of the square, cube ...
% options, level (optional): for the wavelets
% level is an array containing the level of discretization along each
% direction.
%
% Developper : Nicolas Chauffert (nicolas.chauffert@gmail.com)

if nargin==1
    options.h= compute_wavelet_filter('Symmlet',10);
    options.ti=1;
    level=3;
end

d=size(N,2);
if (d==1)
    [pi L]=pi_opt1D(N,options,level);
else
[N_uniq tmp index]=unique(N);
for i=1:size(N_uniq,2)
    [pi L]=pi_opt1D(N_uniq(i),options,level(find(index==i,1)));
    for k = find(index==i)'
       Pi{k}=pi*L;    
    end
end

D=ones(N);
for i=1:d
    dim=ones(1,d);
    dim2=N;
    dim2(i)=1;
    dim(i)=N(i);    
    D=D.*repmat(reshape(Pi{i},dim),dim2);
end
L=sum(D(:));
pi=D/L;
end
end


function [pi L]=pi_opt1D(n,options,level)
% function [pi L]=pi_opt1D(n,options,level)
%
% n=size of the vector
% options, level (optional): for the wavelets
%
% Developper : Nicolas Chauffert (nicolas.chauffert@gmail.com)

Pi=zeros(n,1);
vect=zeros(n,1);
J=log2(n);
for i=1:n
    vect(i)=1;
    FW=perform_wavortho_transfAnis(sqrt(n)*ifft(vect),level,1,options);
    vect(i)=0;
    Pi(i)=max(abs(FW))^2;
end
Pi=fftshift(Pi);
L=sum(Pi(:));
pi=Pi/L;
end
