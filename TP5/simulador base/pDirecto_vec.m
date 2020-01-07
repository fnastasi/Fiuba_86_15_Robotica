function [ pos,config ] = pDirecto_vec(theta_vec,DH)
    
    % largo del vector theta
    n = size(theta_vec,2);
    
    % Inicializo vecor de posición y configuración
    pos = zeros(2,n);
    config = zeros(1,n);
    for i =1:n
       [A_res,config_res] = pDirecto(theta_vec(:,i),DH);
       pos(1,i) = A_res(1,4);
       pos(2,i) = A_res(2,4);
       config(i) = config_res;
    end
    
    
    
    
    
    
    %{
    % Versión anterior
    
    % lango del vector theta
    n = size(theta_vec,2);
    
    % Inicializo vecor de posición y configuración
    pos = zeros(2,n);
    config = zeros(1,n);
    for i =1:n
       t1 =theta_vec(1,i);
       t2 = theta_vec(2,i);
       pos(1,i) = a1*cos(t1) + a2*cos(t1+t2);
       pos(2,i) = a1*sin(t1) + a2*sin(t1+t2);
       config(i) = sign(sin(t2));
    end
    %}
    
    
end