function mask=generateQuadraticSpiralScheme(n,ratio,r0,r1,mask)
%function mask=generateQuadraticSpiralScheme(n,ratio,r0,r1,mask)

if nargin<5
    mask=zeros(n);
end

p=n(1)*n(2)*ratio;
N1=1;
N2=80;
sigma1=max(build_Spiral_quadratic(n,N1,r0,r1),mask);
sigma2=max(build_Spiral_quadratic(n,N2,r0,r1),mask);
while abs(sum(sigma1(:))-sum(sigma2(:)))>50
    N3=N1/2+N2/2;
    sigma3=max(build_Spiral_quadratic(n,N3,r0,r1),mask);
    if sum(sigma3(:))>p
        N2=N3;
    else
        N1=N3;
    end
sigma1=max(build_Spiral_quadratic(n,N1,r0,r1),mask);
sigma2=max(build_Spiral_quadratic(n,N2,r0,r1),mask);
end
mask=sigma2;
end
        
    
function sigma=build_Spiral_quadratic(n,N,r0,r1)
%
%
beta=r0*r1/(r1-r0);
alpha=1;
r=@(x) 1/(1/r0-x/(alpha*beta));
sigma=build_spiral(r,N,n);
end


function sigma= build_spiral(r,N,n)
% r: [0,1] -> R^2
% N an integer (number of revolutions)
% build the spiral s(t)=r(t/N)(cos(2 pi t) , sin(2 pi t))  
%
% Developper : Nicolas Chauffert (2014)

dTheta=1/(2*pi*max(n));
Theta=0;
sigma=zeros(n);
sigma(n(1)/2,n(2)/2)=1;

while (Theta<N)
    Theta=Theta+dTheta;
    s=.5*sqrt(n(1).^2+n(2).^2)*r(Theta/N)*[cos(2*pi*Theta) sin(2*pi*Theta)];
    if round(n(1)/2+s(1))>0 && round(n(1)/2+s(1)) <=n(1) && round(n(2)/2+s(2)) >0 && round(n(2)/2+s(2))<=n(2)
    sigma(round(n(1)/2+s(1)),round(n(2)/2+s(2)))=1;
    end
end
end