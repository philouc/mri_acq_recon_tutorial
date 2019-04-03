function sigma=link2D(p,q,sigma)
p=min(p,size(sigma));
p=max(p,1);
q=min(q,size(sigma));
q=max(q,1);

sigma(q(1),q(2))=1;
sigma(p(1),p(2))=1;

%if(length(p)==2)
if (abs(p(1)-q(1))>abs(p(2)-q(2)))
    if(p(1)-q(1)<0)
        sigma=link2D(q,p,sigma);
    else
        L=p(1)-q(1);
        if (L>1)
           for k=1:(L-1)
               sigma(q(1)+k,round(q(2)+(k/L)*(p(2)-q(2))))=1;
           end               
        end        
    end
else
    if(p(2)-q(2)<0)
        sigma=link2D(q,p,sigma);
    else
        l=p(2)-q(2);
        if (l>1)
           for k=1:(l-1)
               sigma(round(q(1)+(k/l)*(p(1)-q(1))),q(2)+k)=1;
           end               
        end 
    end
end
end