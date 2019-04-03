function sigma = generate2D_CS_scheme(N,ratio,opts)
%
% function sigma = generate2D_CS_scheme(N,ratio,opts)
%
% Creates a 2D sampling scheme of size N:
%
% Input: - N: a 1 x 2 vector of positive integers (size of acquisition
% matrix)
%        - ratio : ratio of measured coefficients
%        - opts .type : sampling strategy ('Independent', 'Radial', 'Radial
%        random', 'Lines', 'Spiral', 'Markov', 'TSP')
%               .distrib : an 2D array of size N, which values are in [0,1]
%               and which sums to 1.
%               .determinist_part: a 2D mask of deterministic measured
%               values.
%
%
% Developper: Nicolas Chauffert, 2014 (nicolas.chauffert@gmail.com)

disp('--------------------')
disp('Generating 2D scheme')
disp('--------------------')

if nargin==1
    ratio=.2;
end
if nargin==2
    opts=struct;
end
    
if ~isfield(opts,'type')
    opts.type='Independent';
end
if ~isfield(opts,'distrib')
    d=distanceMap(N).^(-1);
    opts.distrib=d/sum(d(:));
end   
if ~isfield(opts,'deterministic_part')
    opts.deterministic_part=zeros(N);
    %opts.deterministic_part=set_LF(opts.distrib,2)
end 
opts.distrib(opts.deterministic_part==1)=0;
opts.distrib=opts.distrib/sum(opts.distrib(:));


m=round(N(1)*N(2)*ratio)-sum(opts.deterministic_part(:));
mask=opts.deterministic_part;

if strcmp(opts.type,'Independent')
    opts.distrib=opts.distrib.*(1-opts.deterministic_part);
    opts.distrib=opts.distrib/sum(opts.distrib(:));
    pts=Draw_iid_Points(opts.distrib,m);
    mask(pts(1,:)+(pts(2,:)-1)*N(1))=1;
end

if strcmp(opts.type,'Radial')
    mask=generateSchemeRadial(N,ratio,opts.deterministic_part);
end

if strcmp(opts.type,'Radial Random')
    mask=generateSchemeRadial_RandAngles(N,ratio,opts.deterministic_part);
end

if strcmp(opts.type,'Line')
    l=Draw_iid_Points(sum(opts.distrib,2),round(m/N(2)));
    mask(l,:)=ones(length(l),N(2));
end

if strcmp(opts.type,'Spiral')
    mask=generateQuadraticSpiralScheme(N,ratio,.1,1,opts.deterministic_part);
end

if strcmp(opts.type,'Markov')
    mask=generateSchemeMarkov(opts.deterministic_part, .0001 , m, opts.distrib , 8);
 
end

if strcmp(opts.type,'TSP')
    mask=generateSchemeTravelling(ratio,opts);
end

sigma.mask=mask;
sigma.r=sum(mask(:))/N(1)/N(2);
%figure,imagesc(mask), axis equal, axis off, colormap gray
    