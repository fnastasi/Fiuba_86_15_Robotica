function [M,H,G] = calcMatrizDinamica(X)
global DH;

%Datos de entrada para el problema dinámico
%Parametros cinematicos
a1=DH(1,1);
a2=DH(2,1);
%Parametros dinámicos
global I01zz Xg1 Yg1 m1;
global I02zz Xg2 Yg2 m2;

theta=X(1:2);
theta_p=X(3:4);
g=9.8;

%***************************************************
%CALCULO DE LA MATRIZ M
M(2,2) = I02zz + 2*a2*Xg2*m2 + a2^2*m2;
M(2,1) = M(2,2) + a1* ((a2+Xg2)*cos(theta(2))-Yg2*sin(theta(2)))*m2;
M(1,2) = M(2,1);
M(1,1) = I01zz + 2*a1*Xg1*m1 + a1^2*(m1+m2) + 2*M(1,2) - M(2,2);

% CALCULO DE LA MATRIZ H
H(1,1) = -a1*((a2 + Xg2)*sin(theta(2)) + Yg2*cos(theta(2)))* m2*(2*theta_p(1)*theta_p(2)+theta_p(2)^2);
H(2,1) = -a1*((a2 + Xg2)*sin(theta(2)) + Yg2*cos(theta(2)))* m2*(-theta_p(1)^2);

%CALCULO DE LA MATRIZ G
G(2,1) = m2*g*((Xg2 + a2)*cos(theta(1)+theta(2))- Yg2*sin(theta(1)+theta(2)));
G(1,1) = m1*g*((Xg1 + a1)*cos(theta(1))- Yg1*sin(theta(1))) + m2*g*a1*cos(theta(1)) + G(2,1);
end

