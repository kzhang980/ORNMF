% X = fea;[n,d] = size(X);
% k  =36;
% %E = rand_ortho(n,k);
%  %E = kmeans_ortho(X,k);
% E = sparse(E);
% 
% W0 = E;
% H0 = W0'*X;

%status
% -1: error;
% 0: no potential;
% 1: can reduce
function [W0,e,status,o,T,P] = GOP1(W0,R,l0);


k = size(W0,2);
w = W0(l0,:);
q = find(w>0);
if(length(q)==0); 
    %print('total zero row, warning');
    e = 0;
    status = -1;
    T = 0;
    P = 0;
o = 0;
else
q = q(1);

T = zeros(1,k);
P = zeros(1,k);
for i = 1:k;
    if(i ~=q)
    [obj,~,~,x] = sml(W0,R,q,i,l0);
    else
    [obj,~,x] = sml_0(W0,R,q,l0);
    end;
    T(i) = obj;
    P(i) = x;
end;


dex = find(T<0 & P>0);
if(length(dex) == 0); 
    e = q;
    status = 0;
    o = 0;
else
    mini = min(T(dex));
    e = find(T == mini);
    e = e(1);
    if(e~=q);
        [o,Wpn,Wqn,x] = sml(W0,R,q,e,l0);
        W0(:,e) = Wpn; W0(:,q) = Wqn;
    else
        [o,Wqn,x] = sml_0(W0,R,q,l0);
        W0(:,q) = Wqn;
    end;
    status = 1;
end
%W0 = W0*diag(1./sqrt(1e-10+sum(W0.*W0)));
end;