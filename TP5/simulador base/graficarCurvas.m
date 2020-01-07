function [] = graficarCurvas(time, tita, titaD, titap, titapD, tita2p, tita2pD, u, ptitle)
%graficarCurvas Función para graficar los resultados completos de
%simulaciones, valor medido vs. valor deseado, incluyendo: Torques,
%Posiciones, Velocidades y Aceleraciones, para ambos ejes.
% 
%Parametros:
% time: Vector de tiempo para todos los gráficos
% tita: Matriz de posiciones para ambos ejes, fila 1 corresponde a eje 1
% titaD: Matriz de posiciones deseadas para ambos ejes
% titap: Matriz de velocidades para ambos ejes
% titapD: Matriz de velocidades deseadas para ambos ejes
% tita2p: Matriz de aceleraciones para ambos ejes
% tita2pD: Matriz de aceleraciones deseadas para ambos ejes
% u: Matriz de torques para ambos ejes
% ptitle: Título de la simulación que genera los gráficos

xmax = max(time);

% **************************************************************
% Gráfico de torques vs. posición de ejes
% **************************************************************
plot(time,u);
title(strcat(ptitle, ' - Torques'));
ylabel('Torques [Nm]');
legend({'Eje 1','Eje 2'},'FontSize',7,'Location','northwest');
grid on; xlim([0 xmax]);
xlabel('Tiempo [s]');

% **************************************************************
% Gráfico de posiciones, velocidades y aceleraciones, theta 2
% **************************************************************
for i=1:2
    figure;

    subplot(3,1,1);
    plot(time,(180/pi)*([tita(i,:);titaD(i,:)]));
    title([ptitle ' - Eje ' num2str(i)]);
    ylabel('Pos [°]');
    legend({'Theta','Theta ref'},'FontSize',7,'Location','northwest');
    grid on; xlim([0 xmax]);

    subplot(3,1,2);
    plot(time,[titap(i,:);titapD(i,:)]);
    ylabel('Vel [°/s]');
    legend({'Vel','Vel ref'},'FontSize',7,'Location','northwest');
    grid on; xlim([0 xmax]);

    subplot(3,1,3);
    plot(time,[tita2p(i,:);tita2pD(i,:)]);
    ylabel('Acel [°/s^2]'); 
    legend({'Acel','Acel ref'},'FontSize',7,'Location','northwest');
    grid on; xlim([0 xmax]);
    xlabel('Tiempo [s]');
end

figure;
subplot(2,1,1)
title(ptitle);
plot(time,(u.*titap)');
ylabel('Potencia [W]');
legend({'Eje 1','Eje 2'},'FontSize',7,'Location','northwest');
grid on; xlim([0 xmax]);

subplot(2,1,2)
plot(time,cumsum((u.*titap)'));
ylabel('Energía [J]');
legend({'Eje 1','Eje 2'},'FontSize',7,'Location','northwest');
grid on; xlim([0 xmax]);
xlabel('Tiempo [s]');
end
