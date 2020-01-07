% Modelo en variables de estado para simulación de doble péndulo para TP5
% No incluye el comportamiento elástico de los ejes ni eslabones, por lo
% que los mismos se consideran rígidos  
%
% function dxdt = modeloDinamico(t,x,u)
% 
% t: tiempo actual
% x: vector de estado en el momento actual [theta;thetap] 
% u: variable manipulada, que genera la acción de control

function dxdt = modeloDinamico(t,x,u)

global N;
global Jm;
global Bm;
global Km;

Torq=Km*N*u;
n_ejes=length(Torq);
[M,H,G] = calcMatrizDinamica(x);

thetap = x(n_ejes+1:end);
theta2p = (M+Jm*N^2)\( Torq -Bm*N^2*thetap -H -G);
dxdt = [thetap; theta2p];
end
