% function eval = Linf1_norm()
% 
% This function gathers some function around Linf1 norm
%   Output: eval:
%           eval.function:
%               R^{nd}->R_+: computes the Linf1 norm
%           eval.f_space: computes the norm at each time.
%               R^{nd}->R^n: computes the L1 norm at each time
%           eval.dual: 
%               R^{nd}->R_+: computes the dual norm
%           eval.proxD: 
%               R^{nd}xR_+->R^{nd}: computes the prox of \alpha \|.\|_1inf
%
% Developpers : Pierre Weiss pierre.armand.weiss@gmail.com
%              Nicolas Chauffert nicolas.chauffert@gmail.com

function eval = Linf1_norm()

    
    eval.function= @(x) max(sum(abs(x),2));
    eval.f_space = @(x) sum(abs(x),2);
    eval.dual    = @(x) sum(max(abs(x),[],2));
    eval.proxD   = @(x,alpha) ProxL1_inf(x,alpha);
end

function y=ProxL1_inf(x,alpha)
y=zeros(size(x));
    for i=1:size(x,1)
        y(i,:)=x(i,:)-ProjectionWL1(x(i,:)',ones(size(x,2),1),alpha)';
    end
end
