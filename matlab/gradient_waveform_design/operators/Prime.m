% function sp=Prime(s,Dt)
%
% computes the first derivative of the parametrized curve s.
% In 1D, it corresponds to A*s, where 
%
%        0  0  0  0  ...   0
%       -1  1  0  0  ...   0
%        0 -1  1  0  ...   0
%  A=    .  0 -1  1 0 ..   0
%        .  .  0 -1 1 0.   0
%        .  .  .  ...  ... 0
%        0  0  ...    0 -1 1
%
%
%  Input: s \in R^{nd} (an n x d array)
%
%  Output: sp=A*s \in R^{nd} (an n x d array)
%
% Developpers : Pierre Weiss pierre.armand.weiss@gmail.com
%              Nicolas Chauffert nicolas.chauffert@gmail.com

function sp=Prime(s,Dt)
sp=zeros(size(s));
sp(2:end,:)=s(2:end,:)-s(1:end-1,:);
sp=sp/Dt;