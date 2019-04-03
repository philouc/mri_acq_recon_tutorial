% function eval = Linf2_norm()
% 
% This function gathers some function around L_inf2 norm
%   Output: eval:
%           eval.function:
%               R^{nd}->R_+: computes the L_inf2 norm
%
%               \|s\|_inf2= \| s(i) \|_inf
%
%           eval.f_space: computes the norm at each time.
%               R^{nd}->R^n: computes the L2 norm at each time
%           eval.dual: 
%               R^{nd}->R_+: computes the dual norm
%           eval.proxD: 
%               R^{nd}xR_+->R^{nd}: computes the prox of \alpha \|.\|_12
%
% Inputs are n x d arrays.
%
% Developpers : Pierre Weiss pierre.armand.weiss@gmail.com
%              Nicolas Chauffert nicolas.chauffert@gmail.com

function eval = Linf2_norm()
    
    eval.function=@(x)  max(sqrt(sum(x.^2,2)));
    eval.f_space =@(x)  sqrt(sum(x.^2,2));
    eval.dual=@(x)      sum(sqrt(sum(x.^2,2)));
    eval.proxD=@(x,alpha)proxL12(x,alpha);
end

function y=proxL12(x,alpha)
    eps=1e-10; 
    y=x./repmat(sqrt(sum(x.^2,2)+eps),1,size(x,2)).*max(repmat(sqrt(sum(x.^2,2)),1,size(x,2))-alpha*ones(size(x)),zeros(size(x)));
end

