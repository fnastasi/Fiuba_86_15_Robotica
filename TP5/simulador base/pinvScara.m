function var_art = pinvScara(pose)
    % pose = [x,y,z,g]' ; siendo g el índice de configuración.
    % Las coordenadas de pose deben pasarse en metros
    x = pose(1);
    y = pose(2);
    z = pose(3);
    g = pose(4);
    
    global DH;
    
    % Parámetros del doble péndulo
    a1=DH(1,1); %[m]
    a2=DH(2,1); %[m]
    
    C2 =  (x^2 + y ^2 -a1^2 -a2^2)/(2*a1*a2);
    if abs(C2) <= 1
        S2 = g*sqrt(1-C2^2);
        t2 = atan2(S2,C2);
        
        % Calculo S1 y C1 a partir del sistema de ecuaciones que queda
        P = [x y]';
        M = [a1 + a2*C2, -a2*S2;
            a2*S2, a1+a2*C2];
        res = M\P;
        S1 = res(2);
        C1 = res(1);
        t1 = atan2(S1,C1);
    else
        t2 = NaN;
        t1 = NaN;
    end
    
    var_art = [t1,t2]';
end