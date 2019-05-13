function E = rand_ortho(n,k);
E = zeros(n,k);
 E0 = rand(n,k);
 for i = 1:n;
     [~,idx] = max(E0(i,:));
     E(i,idx) =1;
 end;
  E = E*diag(1./sqrt(sum(E.*E)));