% function sp=Second(s,Dt)
%
% computes -A^T*A*s, where 
%
%        0  0  0  0  ...   0
%       -1  1  0  0  ...   0
%        0 -1  1  0  ...   0
%  A=    .  0 -1  1 0 ..   0
%        .  .  0 -1 1 0.   0
%        .  .  .  ...  ... 0
%        0  0  ...    0 -1 1
%
% is a discrete derivation operator (-A^T*A is a discrete second derivation
% operator).
%
%  Input: s \in R^{nd} (an n x d array)
%
%  Output: sp=A^T*A*s \in R^{nd} (an n x d array)
%
% Developpers : Pierre Weiss pierre.armand.weiss@gmail.com
%              Nicolas Chauffert nicolas.chauffert@gmail.com

function spp=Second(s,Dt)
spp=-PrimeT(Prime(s,Dt),Dt);