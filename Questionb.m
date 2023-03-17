%% Fluid mechanics coursework (Part 2)
% Created by Kenneth Yap
% 21st of November 2021
% Part 2

clear
clc
close all
%% Question a
% Find the optimal d for the maximum N
N = 1:30;
d = 0:0.01:0.5;
count_a = 0;
max_N = 0;
max_d = 0;

% Initial buoyancy
b0 = (273.15 + 17 - 273.15 -7)/(273.15+7)*9.81;
% Quanta emission rate (approximate)
G = 1/60/60; % s^-1

for i = 1:length(N)
    for j = 1:length(d)
    % Constants
    % number of people
    %N = 30; % people

    % depth
    %d = 0.5; %m
    count_a = count_a + 1;

    % ODE for buoyancy
    [t, bcp] = ode45(@func, [0, 3600], [b0 0 0], [], G, d(j), N(i));

    buoyancy = bcp(:,1);
    concentration = bcp(:,2);
    probability = bcp(:,3);

    % Calculation for temperature
    temperature = (buoyancy*(273.15+7)/9.81) + (273.15+7) - 273.15;

    if sum(temperature>=17) == length(temperature) && sum(probability<=0.01) == length(probability)
        max_d = d(j); max_N = N(i);
    end
    end
end

[t, bcp] = ode45(@func, [0, 3600], [b0 0 0], [], G, max_d, max_N);
fprintf('When initial conditions assumed to be zero, d = %g m, the classroom can hold %g students.\n', max_d, max_N)

buoyancy = bcp(:,1);
concentration = bcp(:,2);
probability = bcp(:,3);
temperature = (buoyancy*(273.15+7)/9.81) + (273.15+7) - 273.15;

% Plot graphs
figure
subplot(3,1,1)
plot(t/3600, concentration, 'k')
title('Concentration vs time (zero origin)')
xlabel('time [h]')
ylabel('concentration[m^{-3}]')
ylim([0,0.0015])
xlim([0,1])

subplot(3,1,2)
hold on
plot(t/3600, 0.01*ones(length(t),1),'--r')
plot(t/3600, probability, 'k')
title('Probability vs time (zero origin)')
legend('p = 0.01')
xlabel('time [h]')
ylabel('probability, p')
ylim([0,0.03])
xlim([0,1])
hold off

subplot(3,1,3)
hold on
plot(t/3600, 17*ones(length(t),1),'--r')
plot(t/3600, temperature, 'k')
title('Temperature vs time (zero origin)')
legend('T = 17')
xlabel('time [h]')
ylabel('temperature, T [Celcius]')
ylim([16.9,17.1])
xlim([0,1])
hold off

%% Question b
% Iterate to find the appropriate initial conditions
b0 = (273.15 + 17 - 273.15 -7)/(273.15+7)*9.81;
G = 1/60/60; % s^-1

% Set G to zero when the person leaves the room

count_b = 0;
max_N_b = 0;
max_d_b = 0;
N = 1:30;
d = 0:0.01:0.5;

for i = 1:length(N)
    for j = 1:length(d)
        buoyancy_start = b0;
        concentration_start = 0;
        buoyancy_end = 0;
        concentration_end = 0;
        count_b = count_b + 1;
        % d = 0.5; %m

        for k = 1:100

            % N = 30;
            [~, bcp_iter1] = ode45(@func, [0, 3600], [buoyancy_start concentration_start 0], [], G, d(j), N(i));
            % N = 0;
            [~, bcp_iter2] = ode45(@func, [3600, 3600*1.5], [bcp_iter1(end,1) bcp_iter1(end,2) bcp_iter1(end,3)], [], 0, d(j), 0);

            buoyancy_end = bcp_iter2(end,1);
            concentration_end = bcp_iter2(end,2);

            if (abs(buoyancy_end - buoyancy_start) < 10^-12) && (abs(concentration_end - concentration_start) < 10^-12)
                break
            else
                buoyancy_start = buoyancy_end;
                concentration_start = concentration_end;
            end
        end

        % ODE for buoyancy at N = 30
        %N = 30;
        [t1, bcp_1] = ode45(@func, [0, 3600], [buoyancy_start concentration_start 0], [], G, d(j), N(i));

        buoyancy_1 = bcp_1(:,1);
        concentration_1 = bcp_1(:,2);
        probability_1 = bcp_1(:,3);

        % ODE for buoyancy at N = 0
        %N = 0;
        [t2, bcp_2] = ode45(@func, [t1(end), 3600*1.5], [bcp_1(end,1) bcp_1(end,2) bcp_1(end,3)], [], 0, d(j), 0);

        buoyancy_2 = bcp_2(:,1);
        concentration_2 = bcp_2(:,2);
        probability_2 = bcp_2(:,3);

        % Calculation for temperature
        temperature_1 = (buoyancy_1*(273.15+7)/9.81) + (273.15+7) - 273.15;
        temperature_2 = (buoyancy_2*(273.15+7)/9.81) + (273.15+7) - 273.15;

        % combining array
        t = [t1; t2];
        temperature = [temperature_1; temperature_2];
        concentration = [concentration_1; concentration_2];
        probability = [probability_1; probability_2];
        
    if sum(temperature>=17) == length(temperature) && sum(probability<=0.01) == length(probability)
        max_d_b = d(j); max_N_b = N(i);
    end
    end
end

% ODE for buoyancy at N = 30
% N = 30;

for k = 1:100

    % N = 30;
    [~, bcp_iter1] = ode45(@func, [0, 3600], [buoyancy_start concentration_start 0], [], G, max_d_b, max_N_b);
    % N = 0;
    [~, bcp_iter2] = ode45(@func, [3600, 3600*1.5], [bcp_iter1(end,1) bcp_iter1(end,2) bcp_iter1(end,3)], [], 0, max_d_b, 0);

    buoyancy_end = bcp_iter2(end,1);
    concentration_end = bcp_iter2(end,2);

    if (abs(buoyancy_end - buoyancy_start) < 10^-12) && (abs(concentration_end - concentration_start) < 10^-12)
        break
    else
        buoyancy_start = buoyancy_end;
        concentration_start = concentration_end;
    end
end
        
[t1, bcp_1] = ode45(@func, [0, 3600], [buoyancy_start concentration_start 0], [], G, max_d_b, max_N_b);
fprintf('When initial conditions assumed to be non-zero, d = %g m, the classroom can hold %g students.\n', max_d_b, max_N_b)

buoyancy_1 = bcp_1(:,1);
concentration_1 = bcp_1(:,2);
probability_1 = bcp_1(:,3);

% ODE for buoyancy at N = 0
% N = 0;
[t2, bcp_2] = ode45(@func, [t1(end), 3600*1.5], [bcp_1(end,1) bcp_1(end,2) bcp_1(end,3)], [], 0, max_d_b, 0);

buoyancy_2 = bcp_2(:,1);
concentration_2 = bcp_2(:,2);
probability_2 = bcp_2(:,3);

% Calculation for temperature
temperature_1 = (buoyancy_1*(273.15+7)/9.81) + (273.15+7) - 273.15;
temperature_2 = (buoyancy_2*(273.15+7)/9.81) + (273.15+7) - 273.15;

% combining array
t = [t1; t2];
temperature = [temperature_1; temperature_2];
concentration = [concentration_1; concentration_2];
probability = [probability_1; probability_2];

% Plot graphs
figure

subplot(3,1,1)
hold on
title('Concentration vs time (non-zero origin)')
plot(t/3600, concentration, 'k')
xlabel('time [h]')
ylabel('concentration [m^{-3}]')
ylim([0,0.003])
xlim([0,1.5])
hold off

subplot(3,1,2)
hold on
plot(t/3600, 0.01*ones(length(t),1),'--r')
plot(t/3600, probability, 'k')
legend('p = 0.01')
title('Probability vs time (non-zero origin)')
xlabel('time [h]')
ylabel('probability, p')
ylim([0,0.03])
xlim([0,1.5])
hold off

subplot(3,1,3)
hold on
plot(t/3600, 17*ones(length(t),1),'--r')
plot(t/3600, temperature, 'k')
legend('T = 17')
title('Temperature vs time (non-zero origin)')
xlabel('time [h]')
ylabel('temperature, T [Celcius]')
ylim([12,22])
xlim([0,1.5])
hold off

%% Function to pass into ode45

function final_value = func(t, a, G, d, N)
  b = a(1); 
  c = a(2);
  p = a(3);
  
  % Volume
  V = 200; % m^3
  % heat capacity of air c_p
  c_p = 1.012; % J⋅g−1⋅K−1
  % heat load contribution, Watt
  W = (1000 + N*100); % W
  % buoyancy flux
  F = W*9.81/(c_p*1225*(273.15+7));
  % constant for probability
  r = 0.0002; % m^3/s
  % Constant for vertical opening
  k = 0.5;
  % Width of opening
  w = 6; % m

  Q = k/3*w*sqrt(b.*d.^3);
  dbdt = F/V - Q.*b./V;
  dcdt = G/V - Q.*c./V;
  
  if N>=1
      dpdt = r.*(N-1).*(1-p).*c;
  elseif N == 0
      dpdt = 0;
  end
  
  
  final_value = [dbdt; dcdt; dpdt];
end

