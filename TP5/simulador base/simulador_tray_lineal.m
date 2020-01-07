
% **************************************************************
% Simulador Trayectoria Lineal entre POSE's
% **************************************************************

% Primero se genera un tensor con las matrices homogeneas de las poses

for i =1:size(pdest,2)
    var_art = pinvScara(pdest(:,i));
    [A_aux,config_aux] = pDirecto(var_art,DH);
    POSE(:,:,i) = A_aux;
    POSE_conf(i) = config_aux;
end

% Se comienza calculando la pose inicial. Tener en cuenta que en pdest solo
% se tiene información sobre las coordenadas con respecto a la base de la
% TCP y la configuración. Falta la orientación
%var_art = pinvScara(pdest(:,1));
%[A_aux,config_aux] = pDirecto(var_art,DH);
%POSEA = A; % POSE inicial
%config = config_aux(1);
%phi_actual = 0;

POSEA = POSE(:,:,1); % POSE inicial
%config = POSE_conf(1);
phi_actual = 0;


% Inicializo los estados del sistema
theta = pinvScara(pdest(:,1));
thetap = zeros(size(theta));
theta2p = zeros(size(theta));

pasoSimu = 1;

% Voy a procesar cada una de las instrucciones de movimiento
fprintf('Simulando ');
% Iteración por todos los puntos destinos
for np = 1:size(POSE,3)-1
    config = POSE_conf(np);
    DBA = [ POSE(1:3 ,1:3 ,np)'* POSEA(1:3,1:3) POSE(1:3,1:3,np)'*( POSEA(1:3,4) - POSE(1:3,4,np));
            0 0 0 1];
    DBC = [ POSE(1:3,1:3,np)'* POSE(1:3,1:3,np+1) POSE(1:3,1:3,np)'*( POSE(1:3,4, np+1) - POSE(1:3,4,np));
            0 0 0 1];
    [phi_BA, theta_BA, psi_BA] = invEuler(DBA(1:3,1:3), config, phi_actual);
    [phi_BC, theta_BC, psi_BC] = invEuler(DBC(1:3,1:3), config, phi_BA );
    DA = [DBA(1:3,4) ; theta_BA; psi_BA ];
    DC = [DBC(1:3,4) ; theta_BC; psi_BC ];
    
    % Se calcula tiempo de segmento
    T1 = max([TDeseado(np),2*tacc]);
    
    % Calculo el segmento
   tseg=-tacc+Tm;
   while tseg<=T1-tacc     
      % Indicador de progreso -> 
      if mod(pasoSimu,100)==0
        fprintf('.');
      end
      
      % Obtengo la siguiente referencia usando la función
      % generadorTrayectoriaJoint 
      %[thetaD, thetapD, theta2pD] = interpoladorTrapezoidal(thetaA,thetaB,thetaC,T1,tseg,tacc);
      [THETAD, THETApD, THETA2pD,phi] = interpoladorLineal(DA,DC,phi_BC,phi_BA,T1,tseg,tacc); % THETAD no es el mismo thetaD de la interpolación trapezoidal, es lo de la diapositiva
      D = [dirEuler(phi, THETAD(4,1),THETAD(5,1)), THETAD(1:3,1);0, 0, 0, 1];
      
      %Calculo la POSE en t =tseg
      POS = POSE(:,:,np)*D;
      
      % Calculo el theta deseado con la función inversa manteniendo el
      % índice de configuración
      thetaD = pinvScara([POS(1:3,4); config]);
      
      if pasoSimu > 1
          thetapD = (thetaD - thetaD_ant)/Tm;
          theta2pD = (thetapD - thetapD_ant)/Tm;
      else
          thetapD = zeros(size(thetaD));
          theta2pD = zeros(size(thetaD));
      end
      
      
      %Aplico uno de los siguientes controladores
      
      control_PD;
      %control_PD_peso_propio;
      %control_torque_comp;
      
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
   phi_actual= phi_BC ;
   POSEA=POS; 
end

fprintf(' OK\n');