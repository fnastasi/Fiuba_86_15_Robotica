function [phi, theta, psi] = invEuler(R, g, phi_actual)
    tol = 1E-6;
    
    nx = R(1,1); ny = R(2,1); nz = R(3,1);
    sx = R(1,2); sy = R(2,2); sz = R(3,2);
    ax = R(1,3);ay = R(2,3); az = R(3,3);
    
    if abs(ax) > tol || abs(ay) > tol
        phi = rad2deg(atan2(g*ay,g*ax));
    else
        phi = phi_actual;
    end
    
    theta = rad2deg(atan2(ax*cos(deg2rad(phi)) + ay*sin(deg2rad(phi)),az));
    psi = rad2deg(atan2( -sin(deg2rad(phi))*nx + cos(deg2rad(phi))*ny , -sin(deg2rad(phi))*sx + cos(deg2rad(phi))*sy));
end