function [thetaD, thetapD, theta2pD,phi] = interpoladorLineal(DA,DC,phi_BC,phi_BA,T1,tseg,tacc)
    if tseg < tacc && tseg > -tacc   % Zona 1
        theta2pD = ( (DC/T1) + (DA/tacc) ) /(2* tacc );
        thetapD = (DC/T1)*( tseg + tacc ) /(2* tacc ) + (DA/ tacc )*( tseg - tacc )/(2* tacc);
        thetaD = (DC/T1)*( ( tseg +tacc ) ^2 )/(4* tacc ) + (DA/tacc )* ( ( tseg -tacc ) ^2) /(4*tacc );
        phi = ( phi_BC - phi_BA ) /(2* tacc )*( tseg + tacc )+ phi_BA ;
    elseif tseg >= tacc && tseg < T1 -tacc % Zona 2
        theta2pD = zeros(size(DC));
        thetapD = DC/T1;
        thetaD = (DC/T1)*tseg ;
        phi = phi_BC ;
    end
