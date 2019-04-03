function xx=TSP_sort(pts, concordePath)
% function xx=TSP_sort(pts, concordePath)
% 
% Sort the points in pts (a d x n array of coordonates) using the concorde
% solver.
% concordePath is the path of the concorde folder. (default ./concorde/)
%

n=size(pts,2);
d=size(pts,1);

n_max=10000;
if n>n_max
    xx=[TSP_sort(pts(:,1:n_max),concordePath) TSP_sort(pts(:,(1+n_max):end),concordePath)];
else

%disp('---------------------')
%disp('Computes TSP solution')
%disp('---------------------')
crtFolder=pwd;

%disp('--------------------------------')
%disp('writes distance matrix on a file')
%disp('--------------------------------')
%tic
if nargin == 1
    cd concorde/
else
    cd(concordePath)
end

tmpfile = 'tmp.tsp';
fid = fopen(tmpfile,'W');
if fid<0
    error('Unable to open tmp file.');
    return;
end

% write data to file
fprintf(fid, ['NAME : tmp\nTYPE : TSP\n' ...
    'COMMENT : tmp file from wrapper_linkern.\n' ...
    'DIMENSION : ' num2str(n) '\n' ...
    'EDGE_WEIGHT_TYPE : EXPLICIT\n' ...
    'EDGE_WEIGHT_FORMAT : FULL_MATRIX\n' ... 
    'EDGE_WEIGHT_SECTION\n'] );

for i=1:n
    W=zeros(1,n);
    for j=1:d
    W=W+(pts(j,:)-pts(j,i)).^2;
    end
    W=round(100*sqrt(W));
fprintf(fid, '%d ', W(:)');
end
fprintf(fid, '\nEOF\n');
fclose(fid);
%toc
% call linker

%disp('------------------------------------------')
%disp('computes the solution with CONCORDE solver')
%disp('------------------------------------------')
%tic;
!./LINKERN/linkern -o tmp.res tmp.tsp > tmp.out
%toc;
%% Display
%tic;
tmp=textread('tmp.res');
tmp=tmp(2:end,:);
xx=zeros(d,n);
for i=1:d
xx(i,1:n)=pts(i,tmp(:,1)+1);
end
cd(crtFolder)
end