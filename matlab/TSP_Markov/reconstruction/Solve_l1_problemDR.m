% function rec= Solve_l1_problemDR(data,A,At,niter,gamma)
%
% Solve exact reconstruction with Douglas-Rachford algorithm
%
%        min ||x||_1 s.t Ax=y 
%
%
% input : data  : acquired data (a n x p x q array)
%         A     : a matlab function corresponding to left multiplication by A
%         A     : a matlab function corresponding to left multiplication by At
%         nit   : number of iterations
%         gamma : thresholding parameter (default :1)
%
% output : rec : the reconstructed image (a n x p x q image in the image
% domain)



function rec= Solve_l1_problemDR(y,A,At,niter,gamma)

Norm_y=1;
if nargin <5
    gamma=.05; 
    Norm_y=norm(y(:));
    y=y/Norm_y;
end

Prox_l1 = @(x,tau) max(0,1-tau./max(1e-15,abs(x))).*x;
Proj_set = @(x)x +At(y-A(x));

%Parameters of Douglas Rachford 

lambda=1.5;


z=zeros(size(At(y)));
%L1=zeros(1,niter);

for i=1:niter
    progressbar(i,niter)    
    x=Proj_set(z);
    z=z+lambda*(Prox_l1(2*x-z,gamma)-x);
%    L1(i)=sum(abs(x(:)));
end
%figure, plot(L1)
rec=x*Norm_y;