clc,
close all
clear all

startup;

%% This part defines important parameters to run the code
%
d=3;                % we work with 3D images.
load reference3D      % load I, a reference image volume
N=size(I);
K_space=myFFT(I);   % computes a k space simulation

% Specify wavelet options
options.h= compute_wavelet_filter('Symmlet',10);
myWT = @(x) perform_wavortho_transfAnis(x,[3 3 3],1,options);
myIWT= @(x) perform_wavortho_transfAnis(x,[3 3 3],-1,options);

%% Choose sampling scheme parameters:

% Target distribution
opts.distrib=compute_optimal_distrib(N,options,[3 3 3]);     % Computes optimal distribution

% A deterministic part:
opts.deterministic_part=set_LF(opts.distrib,2);

% Sampling ratio 
R=10;

% Choose the type of sampling scheme, could be:
% TSP or paralel lines

%%
%%%%%%%%%%%%%%%%%%%%%%
%%%% TSP sampling %%%% 
%%%%%%%%%%%%%%%%%%%%%%
mask= generateSchemeTravelling3D(1/R,opts);
figure, plot_3slices(mask)

%% Reconstruct the image from downsampled data

DATA=fftshift(mask).*K_space;
tic,
img=MRI_reconstruct(DATA,fftshift(mask),myWT,myIWT);
toc

%% Display the reconstruction result

figure,
colormap gray
subplot(1,3,1), imagesc(sigma.mask(:,:,round(end/2)), axis equal,axis off
title(['mask: ' opts.type])
subplot(1,3,2), imagesc(img(:,:,round(end/2)), axis equal, axis off
title('reconstructed image')
subplot(1,3,3), imagesc(abs(I(:,:,round(end/2)-img(:,:,round(end/2))), axis equal, axis off
title(['difference with original: (psnr= ' num2str(psnr(img,I)) ')'])

pause()

%%
%%%%%%%%%%%%%%%%%%%%%%
%%% parallel lines %%%
%%%%%%%%%%%%%%%%%%%%%%
opts.distrib=compute_optimal_distrib([N(1) N(2)],options,[3 3]); 
opts.deterministic_part=set_LF(opts.distrib,2);    % Computes optimal distribution
sigma2D=generate2D_CS_scheme(N,1/R,opts);
mask= repmat(sigma2D.mask,[1,1,N(3)]);
figure, imagesc(mask(:,:,round(end/2)))

%% Reconstruct the image from downsampled data

DATA=fftshift(mask).*K_space;
tic,
img=MRI_reconstruct(DATA,fftshift(mask),myWT,myIWT);
toc

%% Display the reconstruction result

figure,
colormap gray
subplot(1,3,1), imagesc(sigma.mask(:,:,round(end/2)), axis equal,axis off
title(['mask: ' opts.type])
subplot(1,3,2), imagesc(img(:,:,round(end/2)), axis equal, axis off
title('reconstructed image')
subplot(1,3,3), imagesc(abs(I(:,:,round(end/2)-img(:,:,round(end/2))), axis equal, axis off
title(['difference with original: (psnr= ' num2str(psnr(img,I)) ')'])
