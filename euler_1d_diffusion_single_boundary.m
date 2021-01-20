L = .1; %length in meters
n = 20; %number of simulation nodes
T0 = 20+273; %initial temperature throughout bar in celcius
t_final = 10000; %simulation time in seconds
dt = .1;% time step in seconds

material = 1; %0=steel, 1=aluminum


%plot selection
animation_plot = false;
over_xl_plot = true;
over_x0_plot = false;
over_time_plot = false;


T1 = 393; %Surface temp 1

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



x = dx/2:dx:L; %Seperate length of bar into mesh elements
t = 0:dt:t_final; %Seperate simulation time into time steps

T = ones(n,1) * T0;
dTdt = zeros(n, 1);

first_node_temp = zeros(2 + (t_final / dt), 1);
end_temp = [T0];
time = 0:(t_final / dt) + 1;

for j = 1:length(t)
    for i = 2:n-1
        dTdt(i) = alpha *(-(T(i) - T(i-1)) / dx^2 + (T(i+1)-T(i)) / dx^2); %Solving for temperature derivitive for interior nodes
    end
    dTdt(1) = alpha *(-(T(1) - T1) / dx^2 + (T(2)-T(1)) / dx^2); %Solving for temperature derivitive for first node
    dTdt(n) = alpha *(-(T(n) - T(n-1)) / dx^2); %Solving for temperature derivitive for final node
    T = T+dTdt*dt;
    newT = T;
    newT(1) = T1;
    first_node_temp(j) = T(1); %Append first node temperature value to list of first node temperature values
    end_temp(end+1) = T(end); %Append end node temperature value to list of end node temperature values
   %time(end+1) = j*dt; 
    
    if animation_plot
        figure(1);
        plot(x, T, 'Linewidth', 3);
        axis([0 L*1.05 min(T1, T0) max(T1, T0)]);
        xlabel(dt* j);
        ylabel('Temperature (\circC)');
    end
end
first_node_temp(1) = T0;
if over_time_plot
    figure(2);
    plot(x, newT, 'Linewidth', 3);
    title(strcat('Temperature across bar at t=', string(t_final), 's'));
    axis([0 L*1.05 -10 100]);
    xlabel('x location (meters)');
    ylabel('Temperature (\circC)');
end
if over_xl_plot
    figure(3);
    plot(time, end_temp, 'Linewidth', 3);
    axis([0, t_final, T0, max(end_temp)]);
    title('Temperature at x=L Over Time')
    xlabel('Time (s)');
    ylabel('Temperature (\circC)');
end
if over_x0_plot
    figure(4);
    plot(time, first_node_temp, 'Linewidth', 3);
    axis([0, t_final, T0, max(first_node_temp)]);
    title('Temperature at x=0 Over Time')
    xlabel('Time (s)');
    ylabel('Temperature (\circC)');
end
%disp(first_node_temp(1))
% T = table(first_node_temp);
% writetable(T,'myDataxl.csv');
% disp(max(first_node_temp))
% disp(min(first_node_temp))
