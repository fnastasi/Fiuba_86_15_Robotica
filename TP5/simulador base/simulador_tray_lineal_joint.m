
% **************************************************************
% Simulador Trayectoria Lineal Joint
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

PxyA = pdest(1:2,1);
PxyB = PxyA;


pasoSimu = 1;

% Voy a procesar cada una de las instrucciones de movimiento
fprintf('Simulando ');
% Iteración por todos los puntos destinos
for np=2:size(pdest,2)
   % Calculo los objetivos del nuevo segmento
   thetaC = pinvScara(pdest(:,np));
   PxyC = pdest(1:2,np);
   DC=thetaC-thetaB;

   % Se calcula tiempo de segmento
   T1 = max([max(DC./v_max ),TDeseado(np),2*tacc]);
   
   
   
   
   
   
   % Calculo el segmento
   tseg=-tacc+Tm;
   while tseg<=T1-tacc     
      % Indicador de progreso -> 
      if mod(pasoSimu,100)==0
        fprintf('.');
      end

      % Obtengo la siguiente referencia usando la función
      % generadorTrayectoriaJoint 
      
      [PxyD, PxypD, Pxy2pD] = interpoladorTrapezoidal(PxyA,PxyB,PxyC,T1,tseg,tacc);
      
      thetaD = pinvScara([PxyD;0;pdest(4,np-1)]);
      if pasoSimu > 1
        thetapD = (thetaD - thetaD_ant)/Tm;
        theta2pD = (thetapD - thetapD_ant)/Tm;
      else
        thetapD = zeros(size(thetaD));
        theta2pD = zeros(size(thetaD));
      end
      
      
      %Aplico uno de los siguientes controladores
            
      %control_PD;
      %control_PD_peso_propio;
      control_torque_comp;
      
      %{
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
     
      %}
      
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
      thetaD_ant = thetaD;
      thetapD_ant = thetapD;
      
      % Proximo paso de simulacion
      pasoSimu=pasoSimu+1;            
      % Calculo el proximo tiempo de segmento
      tseg=tseg+Tm;      
      % Almaceno el proximo tiempo real de la simulación
      tr=tr+Tm;     
   end
   
   % Asignaciones de referencias para calcular el próximo segmento
   %thetaA=thetaD;
   PxyA = PxyD;
   PxyB = PxyC;
   thetaB=thetaC;
end
fprintf(' OK\n');