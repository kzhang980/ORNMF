function [dif,Wqn,x]= sml_0(W0,R,q,l);
 %p=q;
 
    
    
    %o = W0(l,q);
    Wq = W0(:,q);
    Rq = R(:,q);
    o = Wq(l);
    Wqn = Wq;
    
    Wq(l) = 0;
    I = find(Wq>0);
    Wcq = Wq(I);
    Rcq = Rq(I);
    b = Rq(l);
    
    J0 = -(Wcq'*Rcq + o*b);
    
    %Wcqp = Wcq/norm(Wcq,'fro');
    Wcqp = Wcq/sqrt(1-o^2);
    a = Wcqp'*Rcq;
    x = b/sqrt(a^2 + b^2);
    y = sqrt(1-x^2);
    J1 = -(a*y+x*b);
    dif = J1 - J0;

    Wqn(l) = x;
    Wqn(I) = Wcqp*y;
    
