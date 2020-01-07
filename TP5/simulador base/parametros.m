% **************************************************************
% Parametros de la carga
% **************************************************************
ml = 1; % [Kg]
%ml = 0; % [Kg]
Iglzz = 0; % [Kg m^2] la carga está ubicada en el origen de la terna 2

% **************************************************************
% Parametros de los motores seleccionados
% **************************************************************


mm1 = 1.7; % [Kg]
R1 = 0.111/2; % [m]
Igm1zz = (mm1/2)*(R1)^2;  % [Kg m^2]

mm2 = 1.7; % [Kg]
R2 = 0.111/2; % [m]
Igm2zz = (mm2/2)*(R2)^2;  % [Kg m^2]


%{
mm1 = 0; % [Kg]
R1 = 0.111/2; % [m]
Igm1zz = 0;  % [Kg m^2]

mm2 = 0; % [Kg]
R2 = 0.111/2; % [m]
Igm2zz = 0;  % [Kg m^2]
%}

% **************************************************************
% Parametros del robot
% **************************************************************
n_ejes = 2;
global DH;
DH = [0.400 0 0 0;
    0.300 0 0 0]; 
a1=DH(1,1);
a2=DH(2,1);

%Parametros dinamicos del eslabon 1 sin el motor 2 conectado
Igl1zz = 0.07;
Xgl1=-a1/2;
Ygl1=0;
ml1 = 5;

%Parametros dinamicos del eslabon 2 sin la carga útil
Igl2zz = 0.015;
Xgl2=-a2/2;
Ygl2=0;
ml2 = 2;

% Parámetros dinámicos para el eslabón 1 
% calculados al origen de la terna del eslabón. 
% Cuando se incluyan los motores,  suponer que el motor 2
% está fijado al eslabón 1.
global I01zz Xg1 Yg1 m1;
I01zz = (Igl1zz + ml1 * (Xgl1^2+Ygl1^2)) + (Igm2zz);
Xg1 = ((Xgl1*ml1) + (0*mm2))/(ml1+mm2);
Yg1 = ((Ygl1*ml1) + (0*mm2))/(ml1+mm2);
m1 = ml1 + mm2;

global I02zz Xg2 Yg2 m2;
% Parámetros dinámicos para el eslabón 2, calculados al origen de la terna 2
% Cuando se incluya el efecto de la carga, considerarla colgada en el origen
% de la terna 2.
I02zz = (Igl2zz + ml2 * (Xgl2^2+Ygl2^2)) ;
Xg2 = ((Xgl2*ml2) + 0*ml)/(ml2+ml); 
Yg2 = ((Ygl2*ml2)+ 0*ml)/(ml2+ml) ; 
m2 = ml2 + ml;

% **************************************************************
% Parametros del motor
% **************************************************************
% VALORES del manual de Kollmorgen
global Jm;
global Bm;
global N;
global Fm;
global Km;
global fe;
global fm;

% Valores originales
%{
Jm = 1E-5*eye(2,2);      
Bm = 0.000076*eye(2,2);         % Nm/(rad/s)   
N = 1*eye(2,2);
Fm = 0*[2.8/100;2.8/100];	% Si deseara considerar su efecto, incluirla en el modelo

Km=0.10*eye(2);      % Nm/A

% Maximos del motor
v_max = 3000*RPM;   % [rad/seg]
Tau_max = 1;       % [Nm]
%}



%{
% Motor U9D-A
wm = 2*pi/Tm;
%we = ;
Jm = 3.95E-5*eye(2,2); % Nm      
Bm = 5.73E-4*eye(2,2);         % Nm/(rad/s)   
N = 100*eye(2,2);
Fm = 0*[2.8/100;2.8/100];	% Si deseara considerar su efecto, incluirla en el modelo

Km=0.048*eye(2);      % Nm/A

% Maximos del motor
v_max = 6000*RPM;   % [rad/seg]
Tau_max = 3.199;       % [Nm]
%}


% Motor U9D-D
wm = 2*pi/Tm;
%we = ;
Jm = 3.95E-5*eye(2,2); % Nm      
Bm = 0.8*3/(pi*1E3)*eye(2,2);         % Nm/(rad/s)   
N = 100*eye(2,2);
Fm = 0*[2.8/100;2.8/100];	% Si deseara considerar su efecto, incluirla en el modelo

Km=0.076*eye(2);      % Nm/A

% Maximos del motor
v_max = 6000*RPM;   % [rad/seg]
Tau_max = 5.134;       % [Nm]

