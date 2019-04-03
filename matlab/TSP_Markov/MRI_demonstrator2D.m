clc,
close all
clear all

startup;

%% This part defines important parameters to run the code
%
d = 2;              % we work with 2D images.
%load reference      % load I, a reference image
%I=imread('BrainPhantom512.png');
I=imread('BrainPhantom256.png');
N=size(I);
I=double(I);
K_space=myFFT(I);   % computes a k space simulation

% Specify wavelet options
options.h= compute_wavelet_filter('Symmlet',10);
myWT=@(x) perform_wavortho_transf(x,4,1,options);
myIWT=@(x) perform_wavortho_transf(x,4,-1,options);

%% Choose sampling scheme parameters:

% Target distribution
opts.distrib=compute_optimal_distrib(N,options,[3 3]);     % Computes optimal distribution
%opts.distrib=distanceMap(N).^(-2);opts.distrib=opts.distrib/sum(opts.distrib(:));

% A deterministic part:
opts.deterministic_part=set_LF(opts.distrib,4);

% Choose the type of sampling scheme, could be:
% 'TSP', 'Independent', 'Radial', 'Radial_Random',
% 'Markov' or 'Spiral'.

%opts.type='TSP';
opts.type='Independent';
opts.type='Radial';
opts.type='Radial Random';
opts.type='Markov';
opts.type='Spiral';

% Sampling ratio 
R=8;

%% Generate sampling scheme

sigma=generate2D_CS_scheme(N,1/R,opts);

%% Reconstruct the image from downsampled data

DATA=fftshift(sigma.mask).*K_space;
tic,
img=MRI_reconstruct(DATA,fftshift(sigma.mask),myWT,myIWT);
toc

%% Display the reconstruction result

figure
colormap gray
subplot(1,3,1), imagesc(sigma.mask), axis equal,axis off
title(['mask: ' opts.type])
subplot(1,3,2), imagesc(img), axis equal, axis off
title('reconstructed image')
subplot(1,3,3), imagesc(abs(I-img)), axis equal, axis off
title(['difference with original: (psnr= ' num2str(psnr(img,I)) ')'])