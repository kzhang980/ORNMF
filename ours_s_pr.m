function [W0] = ours_s(X,k,T,ratio);

[n,d] = size(X);
W0 = rand_ortho(n,k);

 
for t = 1:T;
    t
  H0 = W0'*X;
  
  R = X*H0';
  W0 = W0*diag(1./sqrt(1e-10+sum(W0.*W0)));
  
  U = randperm(n);
    for j = 1:n;
        j;
        if(rand()>ratio)
            continue;
        end
        [W0,t] = GOP1(W0,R,U(j));
        %count = count + 1;objs(count) = norm(X-W0*H0,'fro').^2;s(count)= 0;
    end
   % W0 = W0*diag(1./sqrt(1e-10+sum(W0.*W0)));
    %W0 = W0.^0.9*diag(1./sqrt(1e-10+sum(W0.^0.9)));
    
 

  
        
end;