function sigma= generateSchemeTravelling3D(r,opts)

Pi=opts.distrib;
sigma=opts.deterministic_part;
N=numel(Pi);
Nf=N*r-sum(sigma(:));
beta=.8;
d=3;
pi=Pi.^(d/(d-1));
pi=pi/sum(pi(:));
Sum=sum(pi(:).^(1-1/d));
n=(Nf/(beta*Sum))^(d/(d-1));
n
n=round(n);
pts=Draw_iid_Points(pi,n);
pts=TSP_sort(pts,opts.Concorde_path);

n=size(pts,2);
N=size(Pi);
for i=1:(n-1)   
    pt1=pts(:,i);
    pt2=pts(:,i+1);
    [X Y Z]=bresenham_line3d(pt1,pt2);
    sigma((Z-1)*N(1)*N(2)+(Y-1)*N(1)+X)=1;  
end
end
