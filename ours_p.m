function Ac = ours_p(X,gnd,k,T,ratio,W0);

[n,d] = size(X);
%W0 = rand_ortho(n,k);
Ac = [];
count = 1;
 
for t = 1:T;
    t;
    
  [~,res] = max(W0');res = bestMap(gnd,res);Ac(count) = length(find(gnd == res))/length(gnd);
  H0 = W0'*X;
  ob(t) = norm(X-W0*H0,'fro');
  R = X*H0';
  W0 = BGOP(W0,H0,R,ratio);
  
  count = count + 1;
   
        
end;

figure,plot(ob);
%hold on,plot(objs,'-');

% figure,
% for i =1:6;
%     for j = 1:6;
%         o = (i-1)*6 + j; s = zeros(90,90); s(:)  = W0(:,o); subplot(6,6,o);imshow(s,[]);colormap(jet)
%     end;
% end;