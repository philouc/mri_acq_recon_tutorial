% function [s, R]=Project_Curve_Affine_Constraints(s0,C,C_linear,Algo_param,R)
%
% This function computes the projection of a curve s0 in R^{nxd} onto the set
% S = {s in R^{nxd}, ||A_k s|| <= alpha_k , 1 <= k <= K}.
%
% INPUT :
% - s0 : curve to project (a n x d array)
% - C : a family of kinematic constraints (see example)
% - C : a family of affine constraints (see example)
% - Algo_param : some algorithm parameters.
% - R (optional argument): initial values of dual variables.
%
% OUTPUT :
% - s : projected curve in R^{nxd} (a n x d array).
% - R : dual variables.
%
%
% Developpers : Pierre Weiss pierre.armand.weiss@gmail.com
%              Nicolas Chauffert nicolas.chauffert@gmail.com

function [s, R]=Project_Curve_Affine_Constraints(s0,C,C_linear,Algo_param,R)

Algo_param=check_params(Algo_param);
[n, d]=size(s0);

n_constraints=size(C,2);
if Algo_param.display_results
    CF=zeros(Algo_param.nit,1);
end

if Algo_param.show_progression
    figure(100)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Set algo parameters %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(Algo_param,'L')
    tic,L=compute_Lipschitz_constant(C,n_constraints,n,d);toc % Need to be fixed ?
else
    L=Algo_param.L;
end
tau = 1/L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Compute dist to constraints %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Algo_param.display_results
    d=zeros(n_constraints, Algo_param.nit);
    for i=1:n_constraints
        d(i,1)=max(C(i).function.function(C(i).operator(s0))-C(i).bound,0);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Pre-compute A_i*s0 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
As=zeros([size(s0) n_constraints]);
for i=1:n_constraints
    As(:,:,i)=C(i).operator(s0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Initialize dual variables %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==5
    if size(R)~=0
        if size(R(:),1)==size(s0(:),1)*n_constraints
            Q=R;
        else
            Q=zeros(size(As));
            R=Q;
            t=sprintf('Input dual variables have wrong size \nInitial dual variables are set to zero');
            disp(t)
        end
    else
        Q=zeros(size(As));
        R=Q;
    end
else
    Q=zeros(size(As));
    R=Q;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% ALGORITHM %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:Algo_param.nit
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Computes \sum_{i=1}^K A_i* s*(i) %%%%%%
    %%%%% (needed to compute the gradient) %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if k>1
        ATq_sum=ATq_sum_last;
    else
        ATq_sum=zeros(size(s0));
        for i=1:n_constraints
            ATq_sum=ATq_sum+C(i).operatorT(Q(:,:,i)) ;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Computes Cost Function %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Algo_param.display_results
        if C_linear.n_linear_constraints>0
            z=s0-ATq_sum;
            s_star=z+C_linear.PI*(C_linear.v-C_linear.A*z);
            CF(k)=sum(s_star(:).*ATq_sum(:))+.5*sum(s_star(:)-s0(:));
            for i=1:n_constraints
                CF(k)=CF(k)-C(i).bound * C(i).function.dual(Q(:,:,i));
            end
        else
            for i=1:n_constraints
                CF(k)=CF(k)-C(i).bound * C(i).function.dual(Q(:,:,i));
            end
            CF(k)=CF(k)-.5*norm(ATq_sum(:))^2 + sum(s0(:).*ATq_sum(:));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Computes the new iterate %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R_prev=R;
    z=s0-ATq_sum;
    if C_linear.n_linear_constraints>0
        s_star=z+C_linear.PI*(C_linear.v-C_linear.A*z);
    else
        s_star=z;
    end
    for i=1:n_constraints
        R(:,:,i)=C(i).function.proxD(Q(:,:,i)+tau*(C(i).operator(s_star)),tau*C(i).bound);
    end;
    Q=R+(k-1)/(k+2)*(R-R_prev);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Computes \sum_{i=1}^K A_i* s*(i) %%%%%%
    %%%(needed at next iteration & to display) %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ATq_sum=zeros(size(s0));
    for i=1:n_constraints
        ATq_sum=ATq_sum+C(i).operatorT(Q(:,:,i)) ;
    end
    ATq_sum_last=ATq_sum;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Compute dist to constraints %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Algo_param.display_results
        for i=1:n_constraints
            d(i,k)=max(C(i).function.function(C(i).operator(s_star))-C(i).bound,0);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Shows algorithm progression %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Algo_param.show_progression
        s=s0;
        for i=1:n_constraints
            s=s-C(i).operatorT(Q(:,:,i));
        end
        Constraints=zeros(size(s,1),n_constraints);
        for i=1:n_constraints
            Constraints(:,i)=C(i).function.f_space(C(i).operator(s));            
            p=2;
            f=figure(100);
            name_fig=['Iteration : ' num2str(k) ' / ' num2str(Algo_param.nit)];
            set(f,'Position',[100 100 500*p 400],'Name',name_fig,'NumberTitle','off');
            subplot(1,p,i);
            plot(s0(:,1),s0(:,2),'k.');
            hold on
            plot_colored_curve(s,Constraints(:,i),C(i).bound)
            hold off
            axis equal
            axis off
            title(['Constraint ' num2str(i)], 'FontSize',44)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the ouptut
z=s0-ATq_sum;
if C_linear.n_linear_constraints>0
    s=z+C_linear.PI*(C_linear.v-C_linear.A*z);
else
    s=z;
end
%%%%%%%%

% display results

Col=['r','g','b','c','m','y','k','w'];
if Algo_param.display_results
    figure(1);
    Constraints=zeros(size(s,1),n_constraints);
    for i=1:n_constraints
        Constraints(:,i)=C(i).function.f_space(C(i).operator(s));
        subplot(2,2,i);
        plot(s0(:,1),s0(:,2),'k.');
        hold on
        plot_colored_curve(s,Constraints(:,i),C(i).bound)
        hold off
        axis equal
        axis off
        title(['Constraint ' num2str(i)], 'FontSize',44)
    end
    
    subplot(2,2,3)
    plot(CF,'linewidth',3);
    l=legend('Cost Function');
    set(l,'FontSize',32)
    xlabel('iteration number')
    title('cost function to MAXIMIZE (dual space)')
    
    subplot(2,2,4)
    for i=1:n_constraints
        hold on
        plot(d(i,:),Col(i),'linewidth',3);
    end
    xlabel('iteration number')
    hold off
    if n_constraints==1
        legend('Constraint 1');
    elseif n_constraints==2
        legend('Constraint 1','Constraint 2');
    elseif n_constraints==3
        legend('Constraint 1','Constraint 2','Constraint 3');
    end
    title('distance to constraints')
    
    
    
    disp('----------------------')
    t=sprintf('Constraint verifications:\n');
    disp(t)
    
    for i=1:n_constraints
        Cons=C(i).function.function(C(i).operator(s));
        disp(['Value of constraint '  num2str(i) ': ' num2str(Cons) '   (Bound: ' num2str(C(i).bound) ')' ])
    end
    disp('----------------------')
    disp('----------------------')
end


end
