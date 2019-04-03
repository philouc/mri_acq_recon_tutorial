% function [z,sigma]=ProjectionWL1(xx,L,alpha)
%
% Computes the projection of xx on a weighted l1 ball
% sigma is the threshold parameter.
% The ball is defined by {x, ||L.*x||_1 \leq \alpha}
%
function [z,sigma]=ProjectionWL1(xx,L,alpha)

x=xx(:);
L=L(:);
n=numel(x);

%Check if solution is non trivial
M=sum(abs(L(:).*x(:)));
if M<=alpha
    %disp('No change')
    z=xx;
    sigma=0;
    return;
end
if alpha<=0
    disp('Are you sure of your parameters?')
    z=zeros(size(xx));
    sigma=infty;
    return;
end

%Compute the projection in the other case...
ax=abs(x(:));
w=ax./L;
[y,J]=sort(w);

E1=L(J(n))*ax(J(n));
E2=L(J(n))^2;
Ep=0;
E=E1-y(n)*E2;
%sprintf('Seuil : %f,E : %f',ax(J(n)),E)
i=n;
while( (E<alpha) && (i>1) )
    i=i-1;
    E1=E1+L(J(i))*ax(J(i));
    E2=E2+L(J(i))^2;
    Ep=E;
    E = E1-y(i)*E2;
    %disp(sprintf('Iteration:%i,Seuil : %f,E : %f',i,ax(J(i)),E));
end

if (i>1)
  a=y(i);
  b=y(i+1);
  r=(Ep-E)/(b-a);
  sigma=(alpha-(Ep-r*b))/r;
else
  sigma=(M-alpha)/E2;   
end

%sprintf('E:%f,Ep:%f,sigma:%f',E,Ep,sigma)
%sigma doit etre compris entre y(i) et y(i+1)

z=x(:);
K=find(ax<sigma*L);
if (~isempty(K))
  z(ax<sigma*L)=0;
end

K=find(ax>=sigma*L);
if (~isempty(K))
  z(K)=x(K)-sigma*sign(x(K)).*L(K);
end

%disp('Change')
z=reshape(z,size(xx));