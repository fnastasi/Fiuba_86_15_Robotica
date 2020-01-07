
% **************************************************************
% Controlador PD
% **************************************************************

% Calculo el torque de control. En este caso es un PD con
% compensación por peso propio
u = Kp*(thetaD-theta)-Kd*thetap;

% Saturación de los motores

Tau_m = Km*u;
%u(abs(u)>Tau_max*[1;1]) = diag(sign(u))*Tau_max*[1;1];
if (abs(Tau_m(1)) > Tau_max)
    u(1) = sign(u(1))*Tau_max/(Km(1,1));
end

if (abs(Tau_m(2)) > Tau_max)
    u(2) = sign(u(2))*Tau_max/(Km(2,2));
end

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
