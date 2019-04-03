function L=compute_Lipschitz_constant(C,n_constraints,n,d)
fct=@f;
function y=f(x)
    y=C(1).operatorT(C(1).operator(reshape(x,n,d)));
    for i=2:n_constraints
        y=y+C(i).operatorT(C(i).operator(reshape(x,n,d)));
    end
    y=y(:);
end
eopts.isreal=false;
L=abs(eigs(fct,n*d,1,'lm',eopts));
end
