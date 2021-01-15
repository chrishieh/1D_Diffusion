L = 0.1; %length in meters
n = 10; %number of simulation nodes
T0 = 1; %initial temperature throughout bar in celcius

T1 = 20; %Surface temp 1
T2 = 100; %Surface temp 2

dx = L/n; %node side in meters
alpha = 0.00001172;

t_final = 20; %simulation time in seconds

dt = 0.1;% time step in seconds

x = dx/2:dx:L-dx/2;

T = ones(n,1) * T0;

dTdt = zeros(n, 1);

t = 0:dt:t_final;

for j = 1:length(t)
    for i = 2:n-1
        dTdt(i) = alpha *(-(T(i) - T(i-1)) / dx^2 + (T(i+1)-T(i)) / dx^2);
    end
    dTdt(1) = alpha *(-(T(1) - T1) / dx^2 + (T(2)-T(1)) / dx^2);
    dTdt(n) = alpha *(-(T(n) - T(n-1)) / dx^2 + (T2-T(n)) / dx^2);
    T = T+dTdt*dt;
    figure(1);
    plot(x, T, 'Linewidth', 3);
    axis([0 L -10 100]);
    xlabel('Distance (m)');
    ylabel('Temperature (\circC)');
    %pause(0.05);
end