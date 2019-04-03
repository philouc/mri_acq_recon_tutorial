% function eval = L1_norm()
% 
% This function gathers some function around L1 norm
%   Output: eval:
%           eval.function:
%               R^{nd}->R_+: computes the L1 norm
%           eval.f_space: computes the norm at each time.
%               R^{nd}->R^n: computes the L1 norm at each time
%           eval.dual: 
%               R^{nd}->R_+: computes the dual norm
%           eval.proxD: 
%               R^{nd}xR_+->R^{nd}: computes the prox of \alpha \|.\|_1
%
% Developpers : Pierre Weiss pierre.armand.weiss@gmail.com
%              Nicolas Chauffert nicolas.chauffert@gmail.com

function eval = L1_norm()
    
    eval.function=@(x) sum(abs(x(:)));
    eval.f_space =@(x) sum(abs(x),2);
    eval.dual=@(x) max(abs(x(:)));
    eval.proxD=@(q,alpha) q-ProjectionWL1(q,ones(size(q)),alpha);
end
