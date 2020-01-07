
% **************************************************************
% Controlador PD con compensación por peso propio
% **************************************************************

% Calculo el torque de control. En este caso es un PD con
% compensación por peso propio
[M,H,G] = calcMatrizDinamica([theta;thetap]);
u = Kp*(thetaD-theta)-Kd*thetap + (N.'*Km)^-1*G ;
%u = Kp*(thetaD-theta)+Kd*(thetapD-thetap) +(Km*N)\(Jeff*theta2pD +Bm*N*N*thetapD) ;


% Integro numericamente sobre el modelo completo de la planta      			
% modeloDinamico es una funcion definida en un script de la siguiente
% manera: 
% function dXdt = modeloDinamico(t,X,Torq)
% Donde t es el tiempo
% X es el vector de estados, es decir [theta;thetap]
% Torq es el vector con el torque de entrada
% Devuelve la derivada del vector de estados (dXdt) que es lo que se
% debe integrar 
% Notar que la acción de control se bloquea durante Tm (es como tener
% un retenedor de orden cero o ZOH)
[tode,X]=ode45(@modeloDinamico,[0 Tm],[theta;thetap],odeOptions,u);

% Lectura de encoder y tacogenerador a partir del modelo integrado
theta=X(end, 1:n_ejes)';
thetap=X(end, n_ejes+1:end)';
