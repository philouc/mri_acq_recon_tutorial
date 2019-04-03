function sigma=generateSchemeMarkov(mask_K, alpha, m, Pi , nb_neighbours)
%function sigma=generateSchemeMarkov(mask_K, alpha, m, Pi , nb_neighbours) 
%
% computes first order markov scheme with center (mask_K) fully centered.
% alpha is the proba of jump
% m is the number of measurements
% Pi is stationnary distribution
% nb_neighboors is 4 or 8
%
% Developper: Nicolas Chauffert (2013) 

sigma=mask_K;
k=sum(sigma(:));
N=size(sigma);
Pi=reshape(Pi,N(1)*N(2),1);

% tirer une paire
while (sum(sigma(:))<m)
    traj=Prand(Pi);
    
    crtLength=1;
    jump=0;
    while (sum(sigma(:))<m && jump==0)
        if (rand()<alpha)
            jump=1;
        end
        if (jump==0)
        last=traj(end);
        voisins=Neighboors(last,N,nb_neighbours,mask_K);
        q=ones(length(voisins),1)/length(voisins); %proposal distribution
        new=voisins(Prand(q));
        accept_rate=min(1,Pi(new)/Pi(last));
        if(rand()<accept_rate) %accept
            traj=[traj new];
            y=mod(new-1,N(1))+1;
            x=(new-y)/N(1)+1;
            if(sigma(y,x)==0)
                k=k+1;
            end
            crtLength=crtLength+1;
            sigma(y,x)=1;
        end
        end

    end
end
end

function voisins=Neighboors(x,N,nb_neigboors,mask_K)

% 4 or 8 neighboors

if (nb_neigboors==8)
    voisins=[x-1 x+1 x-N(1) x+N(1) x-N(1)-1 x-N(1)+1 x+N(1)-1 x+N(1)+1]; % 8 voisins
elseif  (nb_neigboors==4)
    voisins=[x-1 x+1 x-N(1) x+N(1)]; % 4 voisins
end

if (x <=N(1))
    voisins(find(voisins==x-N(1)))=[];
    voisins(find(voisins==x-N(1)-1))=[];
    voisins(find(voisins==x-N(1)+1))=[];
end
if (x>N(1)*(N(2)-1))
    voisins(find(voisins==x+N(1)))=[];
    voisins(find(voisins==x+N(1)-1))=[];
    voisins(find(voisins==x+N(1)+1))=[];
end
if(mod(x,N(2))==1)
    voisins(find(voisins==x-1))=[];
    voisins(find(voisins==x+N(1)-1))=[];
    voisins(find(voisins==x-N(1)-1))=[];
end
if(mod(x,N(2))==0)
    voisins(find(voisins==x+1))=[];
    voisins(find(voisins==x+N(1)+1))=[];
    voisins(find(voisins==x-N(1)+1))=[];
end

if nargin==4
    mask=zeros(N);
    mask(voisins)=1;
    mask=mask.*(1-mask_K);
    voisins=find(mask)';
end
end

function ind = Prand(Pi)
%returns random index between 1 and n (n=size(Pi))
%with the law Pi.

% n=max(size(Pi));
% PiS=zeros(1,n);
% PiS(1)=Pi(1);
% for i=2:n
%     PiS(i)=PiS(i-1)+Pi(i);
% end
PiS=cumsum(Pi);
r=rand();
[m ind]=max(PiS>r);
end