function [dif,Wpn,Wqn,x]= sml(X,W0,H0,R,q,p,l);
 
    o = W0(l,q);
    Wp= W0(:,p);
    Rp = R(:,p);
    Ip = find(Wp>0);
    Wcp = Wp(Ip);
    
    Wq = W0(:,q);
    Rq = R(:,q);
    
    j1 = norm(X-W0*H0,'fro')^2;
%j1 = -2*trace(W0'*R);
    a = Wcp'*Rp(Ip);
    b = Rp(l);
    x = b/sqrt(a^2+b^2+1e-10);
    %J1 = -(Wp'*Rp + Wq'*Rq);
    %J1 = -(a + Wq'*Rq);
    
    y = sqrt(1-x^2);
    Wpn = Wp; Wqn = Wq;
    
    Wpn(l) = x;
    Wpn(Ip) = Wcp*y;
    Wqn(l) = 0;
    Wqn = Wqn/sqrt(1-o^2);% norm(W0(:,q))
    W0(:,p) = Wpn; W0(:,q) = Wqn;
     j2 = norm(X-W0*H0,'fro')^2;
     %j2 = -2*trace(W0'*R);
    
    %J2 = -(a*y + x*Rp(l) + Wqn'*Rq);
dif = j2;%;-(a*(y-1) + x*Rp(l) + (Wqn-Wq)'*Rq);


%     %obj = norm(X-W0*H0,'fro')^2;
%     obj= -trace(W0*R');
%     Z = zeros(1,30);;
% for i = 1:30;
%     Z(i)=W0(:,i)'*R(:,i);
% end;
% ob = -sum(Z);
    
%obj = -(W0(:,p)'*R(:,p) + W0(:,q)'*R(:,q));
    %obj = trace(W0*R');

