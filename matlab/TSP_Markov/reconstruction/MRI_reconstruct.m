function img = MRI_reconstruct(DATA,mask,myWT,myIWT)
%
% function img = MRI_reconstruct(DATA,mask,myWT,myIWT)
%
% MRI reconstructor
%

disp('-----------------------')
disp('Reconstructing MR image')
disp('-----------------------')

A=@(x) mask.*(myFFT(myIWT(x)));
At=@(x) myWT(myIFFT(x));

img= abs(myIWT(Solve_l1_problemDR(DATA,A,At,200)));