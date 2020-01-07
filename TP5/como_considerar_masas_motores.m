% Efecto de las masas de los motores en los parámetros dinámicos del robot
% El siguiente es un cálculo aproximado que desestima el efecto de considerar que el rotor del motor en realidad está girando a una velocidad que es composición de su propia velocidad de rotación y la que le produce el movimiento del eslabón donde está montado. La simplificación es válida porque una velocidad es mucho mayor que la otra (recordar que hay una reducción con N grande)

% Luego el problema se separa en 2. Por un lado consideramos el giro del rotor como si estuviera fijo (con las fórmulas que vimos en clase y que aca no reproducimos) y por otro lado sumamos el efecto de una carga fija al eslabón donde está montado y que es producto principamente de la masa del estator y su distribución.

% Masa del motor 2 que está fijo al eslabón 1
mm2 = 1.7;        % Kg

% Momento de inercia. Como el rotor es muy liviano en comparación al estator, consideramos que el aporte al momento de inercia del eslabón 1 es producido solamente por el estator. Se estima el valor según los datos del manual considerándolo un cilindro de densidad constante
Igm2zz = mm2*(.11/2)^2;  % Kg m^2
% Falta un /2 ?
%Igm2zz = (mm2/2)*(.11/2)^2;  % Kg m^2


% Parámetros dinámicos para el eslabón 1 con la carcaza del motor 2
% calculados al origen de la terna del eslabón. 
global I01zz Xg1 Yg1 m1;

% Momento de inercia del eslabón 1 calculado al origen de la terna 1 (que es donde va montado el motor 2)
% ml1 es la masa del eslabón 1. 
I01zz = (Igl1zz + ml1 * (Xgl1^2+Ygl1^2)) + (Igm2zz);

% Al considerar el motor 2 el centro de masa del eslabón 1 se ve modificado
Xg1 = ((Xgl1*ml1) + (0*mm2))/(ml1+mm2);
Yg1 = ((Ygl1*ml1) + (0*mm2))/(ml1+mm2);

% La "nueva" masa del eslabón 1 incluye al motor 2
m1 = ml1 + mm2;

