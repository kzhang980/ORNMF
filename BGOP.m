% Uratio = 0.3;
% [n,d] = size(X);
% E = rand_ortho(n,k);
% % E = kmeans_ortho(X,k);
% W0 = E;
% H0 = W0'*X;
% R = X*H0';

function W1 = blockGOP(W0,H0,R,Uratio);
W0 = W0*diag(1./sqrt(1e-10+sum(W0.*W0)));
k = size(W0,2);
%J0 = norm(X- W0*H0,'fro')^2;

% current non-zero locations for rows
[~,G] = max(W0');

% searae rows in W0 into c-set (fixed) and u-set (update), make sure that for each column, c set has at least more than 1 non-zeros 
idx_c = [];
idx_u = [];
for i = 1:k;
    ids = find(G == i);
    l = length(ids);
    if(l<3);
        continue;
    end
    %loca = W0(ids,i);
    %[~,odr] = sort(loca,'ascend');
    %ids_perturb = ids(odr);
    ids_perturb = ids(randperm(l));
    lu = max(1,floor(l*Uratio));
    idx_c = [idx_c ids_perturb((lu+1):l)];
    idx_u = [idx_u ids_perturb(1:lu)];
end

% perform GOP on rows in Iu
for i = 1:length(idx_u);
    r = idx_u(i);
    [~,t] = GOP1(W0,R,r);
    if(t>0)
        G(r) = t;
    end;
end;


W1 = zeros(size(W0));
for i = 1:k;
    ic = idx_c(find(G(idx_c)==i));
    iu = idx_u(find(G(idx_u)==i));
    Wic = W0(ic,i);
    Wic = Wic/(1e-10+norm(Wic));
    
    Ric = R(ic,i);
    c = Wic'*Ric;
    u = R(iu,i);
    nu = norm(u,'fro');
    x = u/nu;
    p = c/sqrt(c^2+nu^2+1e-10);
    q = sqrt(1-p^2);
    W1(ic,i) = p*Wic;
    W1(iu,i) = x*q;
end;
W1 = W1*diag(1./sqrt(1e-10+sum(W1.*W1)));
%J1 = norm(X- W1*H0,'fro')^2;