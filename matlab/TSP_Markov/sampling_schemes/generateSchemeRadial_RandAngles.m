function [sigma r]=generateSchemeRadial_RandAngles(N,rt,sigma)
% function [sigma r]=generateSchemeRadial_RandAngles(N,rt,sigma)
%
% Generates a sampling scheme with random angles
% 
center=[N(1)/2+.5 N(2)/2+.5];
if nargin==2
sigma=zeros(N);
end
r=0;

while (r<rt)
    theta=pi*rand()-pi/2;
    if (theta<atan(N(2)/N(1)) && theta >-atan(N(2)/N(1)))
        pt1=[N(1)/2 N(1)/2*tan(theta) ];
        pt2=[-N(1)/2+1 -N(1)/2*tan(theta)+1];
    else
        pt1=[N(2)/2/tan(theta),N(2)/2];
        pt2=[-N(2)/2/tan(theta)+1,-N(2)/2+1];
    end
    sigma=link2D(floor(pt1+center),floor(pt2+center),sigma);
r=sum(sigma(:))/N(1)/N(2);
end
