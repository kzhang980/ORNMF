% %X = Data;
% Chris Ding's method
function [Ac] = nmf0(X,gnd,T);
k = length(unique(gnd));
 [n,d] = size(X);
 W0 = rand(n,k);
 H0 = rand(k,d);
 Ac = [];
 count = 1;
 
for t = 1:T;
    %obj(t) = norm(X - W0*H0,'fro');
    [~,res] = max(W0');
    res = bestMap(gnd,res);
    Ac(count) = length(find(gnd == res))/length(gnd);
    XH= X*H0'; 
    W0 = W0.*XH./(1e-9+W0*(W0'*XH));
    H0 = H0.*(W0'*X)./((W0'*W0)*H0);
    count = count + 1;   
    %S(t) = norm(eye(k) - W0'*W0);
end;

%figure,plot(S);



