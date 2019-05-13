function [Ac,objs] = ours_s(X,gnd,k,T,ratio,W0);

[n,d] = size(X);
%W0 = rand_ortho(n,k);
Ac = [];
count = 1;
 
for t = 1:T;
    t;
  [~,res] = max(W0');res = bestMap(gnd,res);Ac(count) = length(find(gnd == res))/length(gnd);
  H0 = W0'*X;
  count = count + 1;
   objs(count-1) = norm(X-W0*H0,'fro').^2;
  
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
%hold on,plot(objs,'-');

% figure,
% for i =1:6;
%     for j = 1:6;
%         o = (i-1)*6 + j; s = zeros(90,90); s(:)  = W0(:,o); subplot(6,6,o);imshow(s,[]);colormap(jet)
%     end;
% end;