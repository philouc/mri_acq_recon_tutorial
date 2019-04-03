function D = distanceMap(N)
% function D = distanceMap(N)
%
% computes the distance to the center of a discrete grid of size N =[n1 n2
% ...].
%
% Developper : Nicolas Chauffert (2014)
%             nicolas.chauffert@gmail.com
d=length(N);
D=zeros(N);
for i=1:d
    m=((-N(i)/2+.5):(N(i)/2-.5))/N(i);
    d1D=m.^2;
    d_rep=N;
    d_rep(i)=1;
    d_resh=ones(1,d);
    d_resh(i)=N(i);
    D=D+repmat(reshape(d1D,d_resh),d_rep);
end
D=sqrt(D);
end