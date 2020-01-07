function [ A,config ] = pDirecto(theta,DH)
    a1 = DH(1,1);
    a2 = DH(2,1);
    
    t1 = theta(1,1);
    t2 = theta(2,1);
    
    A = [cos(t1+t2),    -sin(t1+t2),    0,  a1*cos(t1) + a2*cos(t1+t2);
         sin(t1+t2),    cos(t1+t2),     0,  a1*sin(t1) + a2*sin(t1+t2);
         0,             0,              1,  0;
         0,             0,              0   1];
    config = sign(sin(t2));
end