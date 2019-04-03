function sigma=generateSchemeRadial(N,r,sigma)
% function sigma=generateSchemeRadial(N,r,sigma)
% sample lines with regular angles
%
%
% Developper: Nicolas Chauffert (nicolas.chauffert@gmail.com), 2014

p1=initial_upper_value(N,r,sigma);
p0=1;
[sigma0 r0]=generateSchemeRadial1(N,p0,sigma);
[sigma1 r1]=generateSchemeRadial1(N,p1,sigma);
while (r0<r && r1>r && (p1-p0)>1)
    p_n=round((p0+p1)/2);
    [sigman r1n]=generateSchemeRadial1(N,p_n,sigma);
    if r1n>r
        p1=p_n;
    else
        p0=p_n;
    end
end
   [sigma r]=generateSchemeRadial1(N,p_n,sigma); 
end


function [sigma r]=generateSchemeRadial1(N,p,sigma)
if (N(2)>N(1))
    [sigma r]=generateSchemeRadial1([N(2) N(1)],p,sigma);
    sigma=sigma';
else
center=[N(1)/2+.5 N(2)/2+.5];
Theta=-pi/2+(1:p)*(pi/p);
%sigma=zeros(N);
for theta=Theta
    if (theta<atan(N(2)/N(1)) && theta >-atan(N(2)/N(1)))
        pt1=[N(1)/2 N(1)/2*tan(theta) ];
        pt2=[-N(1)/2+1 -N(1)/2*tan(theta)+1];
    else
        pt1=[N(2)/2/tan(theta),N(2)/2];
        pt2=[-N(2)/2/tan(theta)+1,-N(2)/2+1];
    end
    sigma=link2D(floor(pt1+center),floor(pt2+center),sigma);
end
r=sum(sum(sigma))/N(1)/N(2);
end
end

function p=initial_upper_value(N,r,sigma)
p=40;
[sigma rc]=generateSchemeRadial1(N,p,sigma);
while rc<r
    p=p+10;
    [sigma rc]=generateSchemeRadial1(N,p,sigma);
end
end
