%function sigma= set_LF(pi,per,connex)
%
% sets a mask of low frequencies, according to the 'per' % of 'pi''s
% largest coefficients; connex is one if one require only a connex area
% (default 1)
%
% Nicolas Chauffert (2014)

function sigma=set_LF(pi,per,connex)
d=nb_dims(pi);
if d==2
[n1 n2 ]=size(pi);
elseif d==3
[n1 n2 n3]=size(pi);
end
piS=sort(pi(:),'descend');
threshold=piS(round(length(piS)*per/100));
sigma=pi>threshold;
 if nargin==2 
     connex=1;
 end
 if connex
     s=regionprops(sigma, 'PixelList');
     i=1;
     while (i<=size(s,1))
         if d==3
         if ismember([n2/2 n1/2 n3/2],[s(i).PixelList(:,1) s(i).PixelList(:,2) s(i).PixelList(:,3)],'rows')
             ind=i;
             i=size(s,1);
         end
         i=i+1;
         elseif d==2
         if ismember([n2/2 n1/2],[s(i).PixelList(:,1) s(i).PixelList(:,2)],'rows')
             ind=i;
             i=size(s,1);
         end
         i=i+1;
         end
     end
     sigma=zeros(size(sigma));
     if d==2
         sigma(min(s(ind).PixelList(:,2)):max(s(ind).PixelList(:,2)),min(s(ind).PixelList(:,1)):max(s(ind).PixelList(:,1)))=1;
     elseif d==3
         sigma(min(s(ind).PixelList(:,2)):max(s(ind).PixelList(:,2)),min(s(ind).PixelList(:,1)):max(s(ind).PixelList(:,1)),min(s(ind).PixelList(:,3)):max(s(ind).PixelList(:,3)))=1;
     end

  if  sum(sigma(:))/numel(pi)*100>per
     while sum(sigma(:))/numel(pi)*100>per
         sigma_last=sigma;
         sigma=imerode(sigma,ones(3*ones(1,d)));
     end
     sigma=sigma_last;
  else
     while sum(sigma(:))/numel(pi)*100<per
         sigma_last=sigma;
         sigma=imdilate(sigma,ones(3*ones(1,d)));
     end
     sigma=sigma_last;
      
 end
end