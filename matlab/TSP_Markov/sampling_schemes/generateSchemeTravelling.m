function sigma= generateSchemeTravelling(r,opts)

Pi=opts.distrib;
sigma=opts.deterministic_part;
N=numel(Pi);
Nf=N*r-sum(sigma(:));
beta=.7;
d=2;
pi=Pi.^(d/(d-1));
pi=pi/sum(pi(:));
Sum=sum(pi(:).^(1-1/d));
n=(Nf/(beta*Sum))^(d/(d-1));
n=round(n);
pts=Draw_iid_Points(pi,n);
pts=TSP_sort(pts,opts.Concorde_path);

%
% To remove:
n=size(pts,2);

for i=1:(n-1)
    pt1=pts(:,i);
    pt2=pts(:,i+1);
    sigma=link2D(pt1',pt2',sigma);    
end
end
