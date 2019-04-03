% function sp=PrimeT(s,Dt)
%
% computes A^T*s, where 
%
%        0  0  0  0  ...   0
%       -1  1  0  0  ...   0
%        0 -1  1  0  ...   0
%  A=    .  0 -1  1 0 ..   0
%        .  .  0 -1 1 0.   0
%        .  .  .  ...  ... 0
%        0  0  ...    0 -1 1
%
% is a discrete derivation operator
%
%  Input: s \in R^{nd} (an n x d array)
%
%  Output: sp=A^T*s \in R^{nd} (an n x d array)
%
% Developpers : Pierre Weiss pierre.armand.weiss@gmail.com
%              Nicolas Chauffert nicolas.chauffert@gmail.com

function sp=PrimeT(s,Dt)
sp=zeros(size(s));
sp(1,:)=-s(2,:);
sp(2:end-1,:)=-s(3:end,:)+s(2:end-1,:);
sp(end,:)=s(end,:);
sp=sp/Dt;
