% function sb=parameterize_maximum_speed(x,v_max,Dt)
%
% computes a parameterized piecewise linear curve at constant speed v_max 
% interpolating the points in list x of size d x n.
%
%
function sb=parameterize_maximum_speed(x,v_max,Dt)
    [d n]=size(x);
    sb=x(:,1);
    r=0;
    d_max=v_max*Dt;
    for i=1:n-1
        crt_vect=x(:,i+1)-x(:,i);
        crt_dist=norm(crt_vect);
        u=crt_vect/crt_dist;
        n_step=floor((crt_dist-r)/d_max);
        if n_step>0
            sb(:,(end+1):(end+n_step+1))=repmat(x(:,i)+r*u,[1 (n_step+1)])+repmat(d_max*u,[1 (n_step+1)]).*repmat(0:(n_step),[d 1]);
            r=d_max-(norm(sb(:,end)-x(:,i+1)));
        else
            r=r-crt_dist;
        end
    end
end