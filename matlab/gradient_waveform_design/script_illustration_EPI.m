%% Gradient waveform design with projection algorithm: the example of EPI
% This document shows how to use the algorithm of trajectory projection to
% design feasible gradient waveforms. It corresponds to Fig. ?? of the paper
% [Chauffert et al., Gradient Waveform Design for variable density sampling
% in Magnetic Resonance Imaging]

%% 
close all
clear all
clc
addpath('MRI/')
addpath('functions/')
addpath('operators/')
folder=strcat(pwd,'/Save_Folder/EPI/Proj80/');
mkdir(folder)

%% Enter the Gradient constraints

% Parameters of the scanner (here use in [Lustig et al, IEEE TMI 2008])
Gmax = 40e-3;  % T/m
Smax = 150e-3; % T/m/ms
Kmax = 500;     % m^-1   

gamma = 42.576*1e3; % kHz/T

alpha = gamma*Gmax;  % in m^-1 ms^-1
beta  = gamma*Smax;  % in m^-1 ms^-2
Dt    = .004;        % sampling time in ms


%% Choose an input trajectory for the algorithm
% Give an EPI trajectory
x=[0;0];
N_line=11;
for i=1:2:N_line
    h1=Kmax*(2*(N_line-i)/(N_line-1)-1);
    h2=Kmax*(2*(N_line-(i+1))/(N_line-1)-1);
    x=[x [[-Kmax Kmax Kmax -Kmax];[h1 h1 h2 h2]]];
end
I=find(abs(x(2,:))>Kmax);
x(:,I)=[];
x=[x [0;0]];
s0=parameterize_maximum_speed(x,.8*alpha,Dt)';
sub=1; % subsampling of the curve for visualization
figure, plot(s0(1:sub:end,1),s0(1:sub:end,2),'b.','linewidth',2)
axis equal, axis off
set(gcf,'Color',[1 1 1])
legend('input trajectory')
saveas(gcf,strcat(folder,'input_traj.pdf'))


[g, param]=exactTSP_pilot(x,Gmax,Smax,gamma,Dt);
T=size(param,1)*Dt;

%% Specify constraints 

%CRV=set_MRI_constraints_RV(alpha,beta,Dt);
CRIV=set_MRI_constraints_RIV(alpha,beta,Dt);

% No additional affine constraint:
C_linear=set_Linear_constraints(size(s0,1),size(s0,2),'start_point',[0 0],'end_point',[0 0]);

Algo_param.nit = 100000;  % number of iterations
Algo_param.L=16;    % Lipschitz constant of the gradient
Algo_param.discretization_step=Dt;

% optional paramteres
Algo_param.show_progression = 0; 
Algo_param.display_results = 0;

%% Project curve with Rotation-Invariant Constraints
tic
s1=Project_Curve_Affine_Constraints(s0,CRIV,C_linear,Algo_param);
toc
% Compute gradients
T_proj=size(s1,1);
t=1:T_proj;
g1=Prime(s1,Dt)/gamma;
gg1=Second(s1,Dt)/gamma;



%% Display projected curves

figure, plot(s0(1:sub:end,1),s0(1:sub:end,2),'--','color',[.5 .5 .5],'linewidth',3)
axis equal, axis off
set(gcf,'Color',[1 1 1])
hold on,
plot(s1(:,1),s1(:,2),'k','linewidth',3)
axis equal, axis off
set(gcf,'Color',[1 1 1])
hold off
title('Projection with RIV constraints')
h=legend('input trajectory', 'projected trajectory','location','sw');
set(h,'FontSize',20,'interpreter','latex');
saveas(gcf,strcat(folder,'proj_traj.pdf'))

%% Display gradients for rotation invariant constraints

T_epi=1.1*T/Dt;
t_epi=1:T_epi;
figure, plot(g(:,1),'r','linewidth',3), axis([0 T_epi -Gmax*1.1 Gmax*1.1])
hold on, plot(g(:,2),'b','linewidth',3), axis([0 T_epi -Gmax*1.1 Gmax*1.1])
%g2n=sqrt(g(:,1).^2+g(:,2).^2);
%hold on, plot(g2n,'k','linewidth',3), axis([0 T_epi -Gmax*1.1 Gmax*1.1])
set(gcf,'Color',[1 1 1])
hold on, plot(t_epi,0*t_epi, '--k','lineWidth',3)
hold on, plot(t_epi,0*t_epi+Gmax, '--k','lineWidth',2)
hold on, plot(t_epi,0*t_epi-Gmax, '--k','lineWidth',2)
set(gca,'XTick',[0,T_epi/1.1])
set(gca,'XTickLabel',{'0',T})
set(gca,'YTick',[-Gmax,Gmax])
set(gca,'YTickLabel',{-Gmax,Gmax})
set(gca,'FontSize',15)
ylabel('T/m')
xlabel('ms')
hold off
h_legend=legend('$ g_x(t)$','$ g_y(t)$');
set(h_legend,'FontSize',20,'interpreter','latex');
title('reparameterization')
saveas(gcf,strcat(folder,'grad_exact.pdf'))

% display gradients
figure, plot(g1(:,1),'r','linewidth',3), axis([0 T_epi -Gmax*1.1 Gmax*1.1])
hold on, plot(g1(:,2),'b','linewidth',3), axis([0 T_epi -Gmax*1.1 Gmax*1.1])
%g1n=sqrt(g1(:,1).^2+g1(:,2).^2);
%hold on, plot(g1n,'k','linewidth',3), axis([0 T_epi -Gmax*1.1 Gmax*1.1])
set(gcf,'Color',[1 1 1])
hold on, plot(t_epi,0*t_epi, '--k','lineWidth',3)
hold on, plot(t_epi,0*t_epi+Gmax, '--k','lineWidth',2)
hold on, plot(t_epi,0*t_epi-Gmax, '--k','lineWidth',2)
set(gca,'XTick',[0,T_proj])
set(gca,'XTickLabel',{'0',T_proj*Dt})
set(gca,'YTick',[-Gmax,Gmax])
set(gca,'YTickLabel',{-Gmax,Gmax})
set(gca,'FontSize',15)
ylabel('T/m')
xlabel('ms')
hold off
h_legend=legend('$ g_x(t)$','$ g_y(t)$');
set(h_legend,'FontSize',20,'interpreter','latex');
title('projection')
saveas(gcf,strcat(folder,'grad_proj.pdf'))

crt_folder=pwd;
cd(folder)
!for a in *.pdf; do pdfcrop "$a"; done
cd(crt_folder)



