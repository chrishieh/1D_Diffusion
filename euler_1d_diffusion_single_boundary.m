L = .1; %length in meters
n = 20; %number of simulation nodes
T0 = 20+273; %initial temperature throughout bar in celcius

material = 1; %0=steel, 1=aluminum

animation_plot = false;
over_xl_plot = true;
over_x0_plot = false;
over_time_plot = false;


T1 = 20+273; %Surface temp 1
T2 = 100+273; %Surface temp 2

dx = L/n; %node side in meters

steel_thermal_conductivity = 50.2; %w/m^2
steel_density = 8050; %kg/m^3
steel_specific_heat = 502; %j/kgk

al_thermal_conductivity = 167;
al_density = 2710;
al_specific_heat = 896;

if material == 0
    alpha = steel_thermal_conductivity / (steel_density * steel_specific_heat);
elseif material == 1
    alpha = al_thermal_conductivity / (al_density * al_specific_heat);
end
disp(alpha);

t_final = 2000; %simulation time in seconds

dt = .1;% time step in seconds

x = dx/2:dx:L+dx/2;

T = ones(n,1) * T0;
size(x)
size(T)
dTdt = zeros(n, 1);

t = 0:dt:t_final;
start_temp = [T0];
end_temp = [T1];
time = [0];
for j = 1:length(t)
    for i = 1:n-1
        if i == 1
            dTdt(i) = alpha * ((T(i+1)-T(i)) / dx^2);
        else
            dTdt(i) = alpha *(-(T(i) - T(i-1)) / dx^2 + (T(i+1)-T(i)) / dx^2);
        end
    end
    %dTdt(1) = dTdt(2);
    dTdt(n) = alpha *(-(T(n) - T(n-1)) / dx^2 + (T2-T(n)) / dx^2);
    T = T+dTdt*dt;
    newT = T;
    newT(n+1) = T2;
    start_temp(end+1) = T(1);
    end_temp(end+1) = T(end);
    time(end+1) = j*dt;
    if animation_plot
        figure(1);
        plot(x, newT, 'Linewidth', 3);
        axis([0 L*1.05 -10 100]);
        xlabel(dt* j);
        ylabel('Temperature (\circC)');
    end
end
if over_time_plot
    figure(2);
    plot(x, newT, 'Linewidth', 3);
    title(strcat('Temperature across bar at t=', string(t_final), 's'));
    axis([0 L*1.05 -10 100]);
    xlabel('x location (meters)');
    ylabel('Temperature (\circC)');
end
if over_x0_plot
    figure(3);
    plot(time, end_temp, 'Linewidth', 3);
    axis([0, t_final, T0, max(end_temp)]);
    title('Temperature at x=0 Over Time')
    xlabel('Time (s)');
    ylabel('Temperature (\circC)');
end
if over_xl_plot
    figure(4);
    plot(time, start_temp, 'Linewidth', 3);
    axis([0, t_final, T0, max(start_temp)]);
    title('Temperature at x=L Over Time')
    xlabel('Time (s)');
    ylabel('Temperature (\circC)');
end
T = table(start_temp);
writetable(T,'C:\Users\Christopher\Desktop\REHTi\1D Diffusion\myDataxl.csv');
disp(max(start_temp))
disp(min(start_temp))
