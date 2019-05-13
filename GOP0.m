function [W0,t] = GOP(W0,H0,R,l0);


k = size(W0,2);
w = W0(l0,:);
q = find(w>0);
if(length(q)==0); 
    %print('total zero row, warning');
    t = 0;
else
q = q(1);

T = zeros(1,k);
P = zeros(1,k);
for i = 1:k;
    if(i ~=q)
    [obj,~,~,x] = sml(W0,H0,R,q,i,l0);
    T(i)=obj;
    P(i) = x;
    end
end;


dex = find(T<0 & P>0);
if(length(dex) == 0); 
    %OB(l+1) = OB(l);
    t = 0;
else
    mini = min(T(dex));
    e = find(T == mini);
    e = e(1);
    [o,Wpn,Wqn,x] = sml(W0,H0,R,q,e,l0);
    W0(:,e) = Wpn; W0(:,q) = Wqn;
    %OB(l+1) = OB(l) + o;    
    t = e;
end
W0 = W0*diag(1./sqrt(1e-10+sum(W0.*W0)));
end;