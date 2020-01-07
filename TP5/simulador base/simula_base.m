%% TP5
%
% Base para resolución TP considerando un control PD
%

clear all
close all
more off

RPM = 2*pi/60;

% **************************************************************
% Instrucciones de movimiento
% **************************************************************
% Matriz de puntos destinos. Los puntos destino se agregan como columnas a la matriz
% Los puntos destinos tienen la forma [x y z g]' donde x,y,z están en m y
% g puede ser 1 o -1
pdest=[-0.300       0.300     0 1;
       0.300        0.300     0 1;
       0.300        0.300     0 1]';
% Los tiempos deseados de cada movimiento se ingresan en un vector columna
TDeseado=[0,    1,  1];

% **************************************************************
% Parametros de la simulacion
% **************************************************************
% Tiempo real
tr=0;
% Tiempo de muestreo
Tm=1E-4;    % s
% Tolerancia del integrador. Parametros de la función ode45
odeOptions = odeset('RelTol',0.001,'AbsTol',0.001,'InitialStep',Tm/10,'MaxStep',Tm/5);

% **************************************************************
% Parametros del robot y de los motores
% **************************************************************
parametros;

% **************************************************************
% Parametros del generador de trayectorias
% **************************************************************
% Tiempo de aceleracion
tacc=100E-3;	% [seg]

% **************************************************************
% Parametros del controlador
% **************************************************************

%Kp=[7000 0;0 2000];
%Kd=[300 0;0 25];

% Valores de la matriz M para realizar aproximación. m22 es constante por
% lo que no se realiza aproximación
t2 = sym('t2');
m22 = I02zz + 2*a2*Xg2*m2 + a2^2*m2;
m12 = m22 + a1* ((a2+Xg2)*cos(t2)-Yg2*sin(t2))*m2; % 6*cos(t2)/25 + 3/20
m12_max = double(vpa(subs(m12, 1)));
m12_min = double(vpa(subs(m12, -1)));
%m12_max = 6/25 + 3/20;
%m12_min = -6/25 + 3/20;
m11_max = I01zz + 2*a1*Xg1*m1 + a1^2*(m1+m2) + 2*m12_min - m22;
m11_min = I01zz + 2*a1*Xg1*m1 + a1^2*(m1+m2) + 2*m12_min - m22;

m11rr = (m11_max + m11_min)/2;
m22rr = I02zz + 2*a2*Xg2*m2 + a2^2*m2;
Jeff = Jm*N*N + [m11rr, 0; 0, m22rr];

wn =2*pi*20 ;
Kp =(Km*N)\(Jeff*wn^2); 
Kp(isnan(Kp)) = 0; % Agrego esta línea porque por algún motivo hacer 0 / "un número" devuelve NaN

Kd = (Km*N)\(2*sqrt(Km*N*Kp*Jeff) - Bm*N*N); 
Kd(isnan(Kd)) = 0; 


simulador_tray_lineal_joint

%simulador_tray_joint;

%simulador_tray_lineal;

%{ 
%Código de back up

% **************************************************************
% Simulador
% **************************************************************
% Siempre arranca desde pdest(:,1) con velocidad nula y termina deteniéndose en el ultimo punto 
% Se usa la función pinvScara, que sirve para calcular el problema inverso
% Calculo los objetivos del nuevo segmento
thetaA = pinvScara(pdest(:,1));
thetaB = thetaA;


% Inicializo los estados del sistema
theta = thetaA;
thetap = zeros(size(theta));
theta2p = zeros(size(theta));


pasoSimu = 1;

% Voy a procesar cada una de las instrucciones de movimiento
fprintf('Simulando ');
% Iteración por todos los puntos destinos
for np=2:size(pdest,2)
   % Calculo los objetivos del nuevo segmento
   thetaC = pinvScara(pdest(:,np));
   DC=thetaC-thetaB;

   % Se calcula tiempo de segmento
   T1 = max([max(DC./v_max ),TDeseado(np-1),2*tacc]);
   
   % Calculo el segmento
   tseg=-tacc+Tm;
   while tseg<=T1-tacc     
      % Indicador de progreso -> 
      if mod(pasoSimu,100)==0
        fprintf('.');
      end

      % Obtengo la siguiente referencia usando la función
      % generadorTrayectoriaJoint 
      [thetaD, thetapD, theta2pD] = interpoladorTrapezoidal(thetaA,thetaB,thetaC,T1,tseg,tacc);
      
      
      % Calculo el torque de control. En este caso es un PD con
      % compensación por peso propio
      u = Kp*(thetaD-theta)-Kd*thetap;
      
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
           
      % Acumulo las variables para poder graficarlas luego
      acum_theta(:,pasoSimu)=theta;
      acum_thetap(:,pasoSimu)=thetap;
      if pasoSimu>1 
          acum_theta2p(:,pasoSimu)=(thetap-thetap_ant)/Tm;% modeloDinamico(0,[theta;thetap],u);
      else
          acum_theta2p(:,pasoSimu)=zeros(size(thetap));
      end
      acum_thetaD(:,pasoSimu)=thetaD;
      acum_thetapD(:,pasoSimu)=thetapD;
      acum_theta2pD(:,pasoSimu)=theta2pD;
      acum_u(:,pasoSimu)=Km*N*u;
      acum_tr(:,pasoSimu)=tr;
      thetap_ant=thetap;
      
      % Proximo paso de simulacion
      pasoSimu=pasoSimu+1;            
      % Calculo el proximo tiempo de segmento
      tseg=tseg+Tm;      
      % Almaceno el proximo tiempo real de la simulación
      tr=tr+Tm;     
   end
   
   % Asignaciones de referencias para calcular el próximo segmento
   thetaA=thetaD;
   thetaB=thetaC;
end
fprintf(' OK\n');
%}


% Calculo las trayectorias
[ pos_ref,config_ref ] = pDirecto_vec(acum_thetaD(:,2:end),DH);
[ pos,config ] = pDirecto_vec(acum_theta(:,1:end-1),DH);

% Grafica de resultados
graficarCurvas(acum_tr, acum_theta, acum_thetaD, acum_thetap, acum_thetapD, acum_theta2p, acum_theta2pD, acum_u, "Pendulo doble");


figure()
plot(pos(1,:),pos(2,:));
hold on
plot(pos_ref(1,:),pos_ref(2,:));
legend('Real','Ref','Location','southwest');
title('Trayectoria');
ylabel('Y [m]');xlabel('X [m]'); grid on
axis('equal');

figure;
delta_pos = pos-pos_ref;
error = sqrt(delta_pos(1,:).^2+delta_pos(2,:).^2);
plot(acum_tr(1:end-1),error*1000);
title('Error de seguimiento');
ylabel('Error [mm]');xlabel('Tiempo [s]'); grid on



figure;
plot(acum_tr(1:end-1),pos_ref(1,:)*1000);
hold on
plot(acum_tr(1:end-1),pos(1,:)*1000);
legend('Ref','Real','Location','southwest');
title('Posición X');
ylabel('Posición X [mm]');xlabel('Tiempo [s]'); 
grid on

figure;
plot(acum_tr(1:end-1),pos_ref(2,:)*1000);
hold on
plot(acum_tr(1:end-1),pos(2,:)*1000);
legend('Ref','Real','Location','southwest');
title('Posición Y');
ylabel('Posición Y [mm]');xlabel('Tiempo [s]'); 
grid on


mean(error)*1000;
std(error)*1000;

if all(eig(N'*Km*Kd + N'*Bm*N) > 0)
    fprintf("La matriz N'KmKd + N'BmN es definida positiva -> Asintoticamente estable\n")
else
    fprintf("La matriz no N'KmKd + N'BmN es definida positiva\n")
end


if any(pos(2,:) > 0.3005) || any(pos(2,:) < 0.2995)
    fprintf('La trayectoria no se encuentra dentro de 0.5 mm de la trayectoria deseada\n');
else
    fprintf('La trayectoria se encuentra dentro de 0.5 mm de la trayectoria deseada\n');
end

