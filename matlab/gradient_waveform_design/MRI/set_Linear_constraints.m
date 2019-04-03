% function C_linear=set_Linear_constraints(n,d,varargin);
%
% Set linear constraints
%
% Input : - n : number of points
%         - d : dimension of signal
%         - (optional) constraints of the form: 'name', value:
%           possible name and values:
%           'start point'             |   [0 0]
%           'end_point'               |   [0 0]
%           'gradient_moment_nulling' |   [1 te] (moment order-echo time)
%           'curve_splitting'         |   100     (in discrete time)
%
% Output : - C_linear: parameters of the affine set of constraints
%
% Example of command:
%
% Developpers  Nicolas Chauffert, Pierre Weiss (2014)

function C_linear=set_Linear_constraints(n,d,varargin)
p=length(varargin)/2;
A=[];
v=[];
n_linear_constraints=0;
if (~isempty(varargin))
    for c=1:2:length(varargin)
        n_linear_constraints=n_linear_constraints+1;
        switch varargin{c}
            case {'start_point'}
                s0=varargin{c+1};
                A=[A ; [1 zeros(1,n-1)]];
                v=[v;s0];
            case {'end_point'}
                se=varargin{c+1};
                A=[A ; [zeros(1,n-1) 1]];
                v=[v;se];
            case {'gradient_moment_nulling'}
                [moment, te]=varargin{c+1};
                T_D=@(x) PrimeT(x',1);
                A=[A ; T_D((1:te).^moment)'];
                v=[v;zeros(1,d)];
            case {'curve_splitting'}
                TR=varargin{c+1};
                k=floor(n/TR);
                M=zeros(k,n);
                for i=1:k
                    M(i,(i-1)*TR+1)=1;
                end
                A=[A;M];
                v=[v;zeros(k,d)];
                
            otherwise
                error(['Invalid optional argument, ', ...
                    varargin{c}]);
        end
    end
    
end

C_linear.n_linear_constraints=n_linear_constraints;
if n_linear_constraints>0
    [M I]=unique(A,'rows');
    A=A(I,:);
    norm=sqrt(sum((A.*A)'))';
    C_linear.A=A./repmat(norm,1,n);
    C_linear.AT=(C_linear.A)';
    C_linear.v=v(I,:)./repmat(norm,1,d);
    C_linear.PI=C_linear.AT*inv(C_linear.A*C_linear.AT);
    C_linear.n_linear_constraints=n_linear_constraints;
end

end