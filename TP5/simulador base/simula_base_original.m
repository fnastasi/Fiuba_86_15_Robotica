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
pdest=[0        -0.6999     1 1;
       0.5576   0.2553      1 1;
       0        -0.6999     1 1;
       0        -0.6999     1 1]';
% Los tiempos deseados de cada movimiento se ingresan en un vector columna
TDeseado=[0 1 0.5 1 1];

% **************************************************************
% Parametros de la simulacion
% **************************************************************
% Tiempo real
tr=0;
% Tiempo de muestreo
Tm=1E-3;    % s
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
Kp=[7000 0;0 2000];
Kd=[300 0;0 25];


% **************************************************************
% Simulador
% **************************************************************
% Siempre arranca desde pdest(:,1) con velocidad nula y termina deteniéndose en el ultimo punto 
% Se usa la función pinvScara, que sirve para calcular el problema inverso
% Calculo los objetivos del nuevo segmento
%thetaA = pinvScara(pdest(:,1));
thetaA = zeros(2,1);
thetaB = thetaA;

% Inicializo los estados del sistema
theta = thetaA;
thetap = zeros(size(theta));

pasoSimu = 1;

% Voy a procesar cada una de las instrucciones de movimiento
fprintf('Simulando ');
% Iteración por todos los puntos destinos
%for np=2:size(pdest,2)
for np=1
    % Calculo los objetivos del nuevo segmento
   %thetaC = pinvScara(pdest(:,np));
   %DC=thetaC-thetaB;

   % Se calcula tiempo de segmento
   %T1 = max([max(DC./v_max ),TDeseado(np),2*tacc]);
   T1=0.5;
   
   % Calculo el segmento
   tseg=-tacc+Tm;
   while tseg<=T1-tacc     
      % Indicador de progreso -> 
      if mod(pasoSimu,100)==0
        fprintf('.');
      end

      % Obtengo la siguiente referencia usando la función
      % generadorTrayectoriaJoint 
      %[thetaD thetapD theta2pD] = interpoladorTrapezoidal(thetaA,thetaB,thetaC,T1,tseg,tacc);
      thetaD=[pi/4 pi/4]'; 
      thetapD=zeros(2,1);
      theta2pD=zeros(2,1);
      
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
   %thetaA=thetaD;
   %thetaB=thetaC;
end
fprintf(' OK\n');

% Calculo las trayectorias
%[ pos_ref,config_ref ] = pDirecto(acum_thetaD(:,2:end),DH);
%[ pos,config ] = pDirecto(acum_theta(:,1:end-1),DH);

% Grafica de resultados
graficarCurvas(acum_tr, acum_theta, acum_thetaD, acum_thetap, acum_thetapD, acum_theta2p, acum_theta2pD, acum_u, "Pendulo doble");


%figure()
%plot(pos(:,1),pos(:,2));
%hold on
%plot(pos_ref(:,1),pos_ref(:,2));
%legend('Ref','Real','Location','southwest');
%title('Trayectoria');
%ylabel('Y [m]');xlabel('X [m]'); grid on
%axis('equal');

%figure;
%delta_pos = pos-pos_ref;
%error = sqrt(delta_pos(:,1).^2+delta_pos(:,2).^2);
%plot(acum_tr(1:end-1),error*1000);
%title('Error de seguimiento');
%ylabel('Error [mm]');xlabel('Tiempo [s]'); grid on

%mean(error)*1000
%std(error)*1000
