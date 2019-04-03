function f = perform_wavortho_transfAnis(f,level,dir,options)

% perform_wavortho_transfAnis - compute orthogonal wavelet transform
%
%   fw = perform_wavortho_transfAnis(f,level,dir,options);
%
%   You can give the filter in options.h.
%
%   Works in arbitrary dimension.
%
%   Copyright (c) 2009 Gabriel Peyre
%
%  modified: 2014 Nicolas Chauffert
%
% level is an array of values : level(d) is the number of levels in the
% dimension d.
%
%



options.null = 0;
h = getoptions(options,'h', compute_wavelet_filter('Daubechies',4) );
g = [0 h(length(h):-1:2)] .* (-1).^(1:length(h));

N=size(f);

for i=1:nb_dims(f)
    if mod(N(i),2^level(i))~=0
        error ('size should be a multiple of level(i)')
    end
end

if dir==1
    %%% FORWARD %%%
    for j=1:max(level)
        for i=1:nb_dims(f)
            W(i) = N(i)/2^min(j-1,level(i));
        end
        a = subselect_d(f,W);
        for d=1:nb_dims(f)
            if j<=level(d)
             a = cat(d, subsampling(cconv(a,h,d),d), subsampling(cconv(a,g,d),d) );
            end
        end
        f = subassign(f,W,a);
    end
else
    %%% BACKWARD %%%
    for j=max(level):-1:1
        for i=1:nb_dims(f)
            W(i) = N(i)/2^min(j-1,level(i));
        end
        a = subselect_d(f,W);
        for d=1:nb_dims(f)
            if j<=level(d)
             w = subselectdim(a,W(d)/2+1:W(d),d);
             a = subselectdim(a,1:W(d)/2,d);
             a = cconv(upsampling(a,d),reverse(h),d) + cconv(upsampling(w,d),reverse(g),d);
            end
        end
        f = subassign(f,W,a);
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = subselect_d(f,W)
dim=size(W,2);
if dim==1
    f = f(1:W(1));
elseif dim==2
    f = f(1:W(1),1:W(2));
elseif dim==3
    f = f(1:W(1),1:W(2),1:W(3));
else
    error('Not implemented');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = subselectdim(f,sel,d)
switch d
    case 1
        f = f(sel,:,:,:,:,:,:,:);
    case 2
        f = f(:,sel,:,:,:,:,:,:);
    case 3
        f = f(:,:,sel,:,:,:,:,:);
    case 4
        f = f(:,:,:,sel,:,:,:,:);
    case 5
        f = f(:,:,:,:,sel,:,:,:);
    case 6
        f = f(:,:,:,:,:,sel,:,:);
    case 7
        f = f(:,:,:,:,:,:,sel,:);
    case 8
        f = f(:,:,:,:,:,:,:,sel);
    otherwise
        error('Not implemented');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = subassign(f,W,g)
dim=size(W,2);
if dim==1
    f(1:W(1)) = g;
elseif dim==2
    f(1:W(1),1:W(2)) = g;
elseif dim==3
    f(1:W(1),1:W(2),1:W(3)) = g;
else
    error('Not implemented');
end
