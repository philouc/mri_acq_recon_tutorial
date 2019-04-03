% function C = set_MRI_constraints_RIV(alpha,beta,Dt)
%
% Set gradient constraints (rotation invariant)
%
% INPUT : - alpha: maximal gradient value (in cm^-1.ms^-1)
%         - beta : gradient slew_rate (in cm^-1.ms^-2)
%         - Dt : discretization step of the cruve
%
% OUTPUT : - C : contraints
%
% Developper: Nicolas Chauffert (nicolas.chauffert@gmail.com)


function C = set_MRI_constraints_RIV(alpha,beta,Dt)
% Constraint 1:
C1.function=Linf2_norm;
C1.operator=@(s) Prime(s,1);
C1.operatorT=@(s) PrimeT(s,1);
C1.bound=alpha*Dt;

% Constraint 2:
C2.function=Linf2_norm;
C2.operator=@(s) Second(s,1);
C2.operatorT=@(s) Second(s,1);
C2.bound=beta*Dt^2;

C=[C1 C2];
