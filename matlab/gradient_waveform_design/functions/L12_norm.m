% function eval = L12_norm()
% 
% This function gathers some function around Linf norm
%   Output: eval:
%           eval.function:
%               R^{nd}->R_+: computes the L12 norm
%           eval.f_space: computes the norm at each time.
%               R^{nd}->R^n: computes the L2 norm at each time
%           eval.dual: 
%               R^{nd}->R_+: computes the dual norm
%           eval.proxD: 
%               R^{nd}xR_+->R^{nd}: computes the prox of \alpha \|.\|_inf2
%
% Developpers : Pierre Weiss pierre.armand.weiss@gmail.com
%              Nicolas Chauffert nicolas.chauffert@gmail.com

function eval = L12_norm()
    
    eval.function=@(x) sum(sqrt(sum(x.^2,2)));
    eval.f_space =@(x) sqrt(sum(x.^2,2));
    eval.dual=@(x) max(sqrt(sum(x.^2,2)));
    eval.proxD=@(q,alpha) proxLinf2(q,alpha);
end

function y=proxLinf2(q,alpha)
    eps=1e-10;
    Q_norm=sqrt(sum(q.^2+eps,2));
    Coeffs=Q_norm-ProjectionWL1(Q_norm,ones(size(Q_norm)),alpha);
    y=repmat(Coeffs./Q_norm,1,size(q,2)).*q;
end
