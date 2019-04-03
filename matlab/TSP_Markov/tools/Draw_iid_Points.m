function pts=Draw_iid_Points(p,n)
% function pts=Draw_iid_Points(p,n)
%
% Inputs :   p: a pdf (array of dimension d)
%            n: number of points to draw independently
%
% Outputs :  pts: list of points drawn

dimension=length(size(p));
if (dimension==2)
    if (size(p,1)==1 || size(p,2)==1)
        dimension=1;
    end
end
distrib=cumsum(p(:));

pp=rand(n,1);
pp=sort(pp(:));

pts_i=zeros(n,1);
i=1;
iprec=1;
ind=1;
while (i<=n)
    while (distrib(ind)<pp(i))
        ind=ind+1;            
    end
    while (i<=n && pp(i)<=distrib(ind))
        i=i+1;
    end
    pts_i(iprec:i-1)=ind;    
    iprec=i;
    ind=ind+1;
   % if ind>numel(distrib)
   %     i=n+1;
   % end
end
%pts_i=unique(pts_i);

pts=zeros(dimension,n);
if dimension==1
    pts=pts_i;
elseif dimension==2
    [p1,p2]=ind2sub(size(p),pts_i);
    pts(1,:)=p1;pts(2,:)=p2;
elseif dimension==3
    [p1,p2,p3]=ind2sub(size(p),pts_i);
    pts(1,:)=p1;pts(2,:)=p2;pts(3,:)=p3;
end

pts=unique(pts','rows')';
if dimension==1
    pts=unique(pts)';
end
% refine step:
if dimension>1
    n_new=n-size(pts,2);
elseif dimension==1
    n_new=n-length(pts);
end

if (n_new>1)
    if dimension==1
        p(pts)=0;
    elseif dimension>1
    p(pts_i)=0;
  
    end
    pts=[pts Draw_iid_Points(p/sum(p(:)),n_new)];
end
