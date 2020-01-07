function [thetaD, thetapD, theta2pD] = interpoladorTrapezoidal(thetaA,thetaB,thetaC,T1,tseg,tacc)
    DA = thetaA - thetaB;
    DC = thetaC - thetaB;
    
    % Zona 1
    if tseg <= tacc && tseg >= -tacc
        thetaD =(DC/T1)*((tseg+tacc)^2)/(4*tacc) + (DA/tacc)*((tseg-tacc)^2)/(4*tacc) +thetaB;
        thetapD = DC*( (tseg+tacc)/(2*tacc))/T1 + DA*((tseg-tacc)/(2*tacc))/tacc; 
        theta2pD =  (DC/T1 + DA/tacc) / (2*tacc);
    % Zona 2
    elseif tseg >= tseg && tseg <= T1 -tacc
        thetaD =(DC/T1)*tseg +thetaB;
        thetapD = DC/T1; 
        theta2pD =  0;
    else
        % Para indicar que se ingreso un valor de tseg que no es coherente.
        thetaD =NaN;
        thetapD = NaN; 
        theta2pD =  NaN; 
    end