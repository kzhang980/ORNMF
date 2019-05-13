function [dif,y1,y2,x]= sml(W0,H0,R,q,p,l);
 
    o = W0(l,q);
    Wp= W0(:,p);
    Rp = R(:,p);
    Ip = find(Wp>0);
    Wcp = Wp(Ip);
    
    Rq = R(:,q);
    Wq = W0(:,q);
    J1 = -(Wp'*Rp + Wq'*Rq);

    a = Wcp'*Rp(Ip);
    b = Rp(l);
    x = b/sqrt(a^2+b^2+1e-10);
    
    W0(l,p) = x;
    W0(Ip,p) = Wcp*sqrt(1-x^2);
    W0(l,q) = 0;
    W0(:,q) = W0(:,q)/sqrt(1-o^2);% norm(W0(:,q))
    y1 = W0(:,p);
    y2 = W0(:,q);
    J2 = -(W0(:,p)'*Rp + W0(:,q)'*Rq);
dif = J2 - J1;




%     %obj = norm(X-W0*H0,'fro')^2;
%     obj= -trace(W0*R');
%     Z = zeros(1,30);;
% for i = 1:30;
%     Z(i)=W0(:,i)'*R(:,i);
% end;
% ob = -sum(Z);
    
%obj = -(W0(:,p)'*R(:,p) + W0(:,q)'*R(:,q));
    %obj = trace(W0*R');

