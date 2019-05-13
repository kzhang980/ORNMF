function W1 =RBGOP(W0,R,ratio);
W0 = W0*diag(1./sqrt(1e-10+sum(W0.*W0)));
[n,k] = size(W0);

% current non-zero locations for rows
[~,G] = max(W0');
G = G';

% searae rows in W0 into c-set (fixed) and u-set (update), make sure that for each column, c set has at least more than 1 non-zeros 
idx_c = [];
idx_u = [];
for i = 1:k;
    ids = find(G == i);
    l = length(ids);
    if(l<3);
        continue;
    end
    loca = W0(ids,i);
    [~,odr] = sort(loca,'ascend');
    ids_perturb = ids(odr);
    ids_perturb = ids(randperm(l));
    lu = max(1,floor(l*ratio));
    idx_c = [idx_c; ids_perturb((lu+1):l)];
    idx_u = [idx_u; ids_perturb(1:lu)];
end
C = zeros(n,1);
U = zeros(n,1);
C(idx_c) = 1;
U(idx_u) = 1;


W1 = zeros(size(W0));
for i = 1:k;
    ic = find(C == 1 & G == i);

    iu = find(U == 1 & G == i);
    
    Wic = W0(ic,i);
    Wic = Wic/(1e-10 +norm(Wic));
    
    Ric = R(ic,i);
    c = Wic'*Ric;
    u = R(iu,i);
    nu = norm(u,'fro');
    x = u/(1e-10+nu);
    p = c/sqrt(1e-10+c^2+nu^2);
    q = sqrt(1-p^2);
    W1(ic,i) = p*Wic;
    W1(iu,i) = x*q;
end;
