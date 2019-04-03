% function eval = Linf_norm()
% 
% This function gathers some function around Linf norm
%   Output: eval:
%           eval.function:
%               R^{nd}->R_+: computes the Linf norm
%           eval.f_space: computes the norm at each time.
%               R^{nd}->R^n: computes the Linf norm at each time
%           eval.dual: 
%               R^{nd}->R_+: computes the dual norm
%           eval.proxD: 
%               R^{nd}xR_+->R^{nd}: computes the prox of \alpha \|.\|_1
%
% Developpers : Pierre Weiss pierre.armand.weiss@gmail.com
%              Nicolas Chauffert nicolas.chauffert@gmail.com

function eval = Linf_norm()
    
    eval.function=@(x) max(abs(x(:)));
    eval.f_space =@(x) max(abs(x),[],2);
    eval.dual=@(x) sum(abs(x(:)));
    eval.proxD=@(q,alpha) sign(q).*max(abs(q)-alpha,0);
end
