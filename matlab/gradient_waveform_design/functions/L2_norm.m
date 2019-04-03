% function eval = L2_norm()
% 
% This function gathers some function around L2 norm
%   Output: eval:
%           eval.function:
%               R^{nd}->R_+: computes the L2 norm
%           eval.f_space: computes the norm at each time.
%               R^{nd}->R^n: computes the L2 norm at each time
%           eval.dual: 
%               R^{nd}->R_+: computes the dual norm
%           eval.proxD: 
%               R^{nd}xR_+->R^{nd}: computes the prox of \alpha \|.\|_2
%
% Developpers : Pierre Weiss pierre.armand.weiss@gmail.com
%              Nicolas Chauffert nicolas.chauffert@gmail.com

function eval = L2_norm()

    eps=1e-10; 
    
    eval.function= @(x) sqrt(sum(x(:).*x(:)));
    eval.f_space = @(x) sqrt(sum(x.*x,2));
    eval.dual    = @(x) sqrt(sum(x(:).*x(:)));
    eval.proxD   = @(q,alpha) max(q.*(1-alpha/sqrt(sum(q(:).*q(:))+eps)),0);
end
