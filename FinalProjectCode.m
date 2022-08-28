%% CASE 600 - Set up
% MAE 600 Spring 2021
% Matthew Nigro 


% Read TMY weather data
data = csvread('DryBulbTemperature(1).csv',1,0);
Lat = 39.76*pi/180;    % degrees
Long = -104.86*pi/180; % degrees

% Set sampling interval for simulation
dt = 5*60; % Sampling interval - 1 or 5 minute

% Building model

% Building dimensions
Asouth = 8*2.7 - 2*2*3;
Anorth = 8*2.7;
Aeast = 6*2.7;
Awest = 6*2.7;
Aroof = 6*8;
Awin = 2*2*3;

% R-values
% R' = 1/UA
Uwall = 32.715;
Uroof = 15.253;
Utot = Uwall + Uroof;
Rtot = 1/Utot;
R1 = Rtot/2;
R2 = R1; 

% % R-values - another approach to calculate R values
% % R'=1/R
% Rwall_unit = 0.121+0.075+1.65+0.064+0.034;
% Rwall = Rwall_unit /  (Asouth+Anorth+Aeast+Awest);
% Rroof_unit = 0.121+0.063+2.794+0.136+0.034;
% Rroof = Rroof_unit / Aroof;
% Rtot = 1/ (1/Rwall + 1/Rroof); 
% R1 = Rtot / 2;

Uwin = 36;
Rwin = 1/Uwin;

% C-values
% C = (rho)(Cp)(Vol)
rho_air = 1.225 * 0.8;  
cp_air = 1005;  % J/(kg*K)
vol_air = 6 * 8 * 2.7;
Cz = rho_air * cp_air * vol_air / 1;

Cwall = (950*840*0.012 + 12*840*0.066 + 530*900*0.009) * (Asouth+Anorth+Aeast+Awest);  % J/K
Croof = (950*840*0.010 + 12*840*0.1118 + 530*900*0.019) * Aroof;
Cw = Cwall + Croof;  % J/K

% HVAC coefficient of performance
cop = 1.0; 

% System matrices
A = [-(1/Cw)*(1/R1+1/R2), 1/(Cw*R2); ...
    1/(Cz*R2), -(1/Cz)*(1/R2+1/Rwin)];
Ad = expm(A*dt); % e^(A*dt)

B = [0; ...
    cop/Cz];
Bd = inv(A)*(Ad-eye(2))*B;

E = [1/(Cw*R1), 1/Cw, 0; ...
    1/(Cz*Rwin), 0, 1/Cz];

Ed = inv(A)*(Ad-eye(2))*E;

% Disturbances
% Interpolate ambient temperature
Tamb_hourly = data;
t_hourly = 1/(3600/dt):length(Tamb_hourly);
t_interp = 1/(3600/dt):dt/60/60:length(Tamb_hourly);
Tamb_interp = interp1(t_hourly,Tamb_hourly,t_interp,'linear','extrap');

% Solar radiation heat gains through 4walls + roof + window
% StartDate = datetime(1990,1,1,00,00,00);
% EndDate = datetime(1991,1,1,00,00,00);
% SolSouth = SolRad(dt,Lat,pi,1,0, StartDate, EndDate) * Asouth;
% solNorth = SolRad(dt,Lat,0,1,0, StartDate, EndDate) * Anorth ;
% SolEast = SolRad(dt,Lat,pi/2,1,0, StartDate, EndDate) * Aeast ;
% SolWest = SolRad(dt,Lat,3*pi/2,1,0, StartDate, EndDate) * Awest;
% SolRoof = SolRad(dt,Lat,0,0,0, StartDate, EndDate) * Aroof;
% SolWindow = SolRad(dt,Lat,pi,1,1, StartDate, EndDate) * Awin;
 
% SolWall = SolSouth + solNorth + SolEast + SolWest + SolRoof + SolWindow;
% save the solar radiation data
% it takes a long time to calculate
%save('case_600_data')
load('case_600_data')

% Constant internal heat gains
Qint = 100*[ones(8*60*60/dt,1); zeros(10*60*60/dt,1); ones(6*60*60/dt,1)]; % out for work during the day from 8-18
% plot(Qint)
Qint = repmat(Qint,365,1);  % repmat([10;11],3,1); repeat matrix
% Collect disturbances
W = [Tamb_interp; SolWall'./100; Qint'];

% OPEN LOOP SIMULATION

% control step
control_dt = 15*60;  % in seconds, SolWall is every dt second
n = 1; 

% Initialize states
x0 = [22; 22];  % initial temperature values
x = zeros(5,length(SolWall)/(control_dt/dt)); % Twall, Tzone and Tamb
x(1:2,1) = x0;

% Simulation w/ zero HVAC power
for t = 1:control_dt/dt:length(SolWall) % 1:length(SolWall)
    x(1:2,n+1) = Ad*x(1:2,n) + Bd*0 + Ed*W(:,t);
    x(3:5,n) = W(:,t);
    n = n + 1;  % n is used to store x data, x and W have different time scale
end

% Plot results
figure; plot(x(1:3, :)');
xticks([24*1*4 24*32*4 24*60*4 24*91*4 24*121*4 ...
    24*152*4 24*182*4 24*213*4 24*244*4 24*274*4 24*305*4 24*335*4 24*365*4])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan'})
xlabel('Month')
ylabel('Temperature (\circ C)')
legend({'Twall','Tzone','Tamb'}); title('Open Loop Simulation')

% CLOSED LOOP SIMULATION - ON/OFF Control

% control step
control_dt = 15*60;  % in seconds, SolWall is every dt second
n = 1; 

% Initialize states
xx = zeros(5,length(SolWall)/(control_dt/dt)); % Twall, Tzone and Tamb
xx(1:2,1) = x0; 


u = zeros(1,length(SolWall)/(control_dt/dt));
pHVAC = 2000;

Tz_hi = [24*ones(8*12,1); 26*ones(10*12,1); 24*ones(6*12+1,1)];
Tz_hi = repmat(Tz_hi,365,1); 

Tz_lo = [20*ones(8*12,1); 18*ones(10*12,1); 20*ones(6*12+1,1)];
Tz_lo = repmat(Tz_lo,365,1); 

%Summer setpoints
Tz_hi_summer = [25*ones(8*12,1); 27*ones(10*12,1); 25*ones(6*12+1,1)];
Tz_hi_summer = repmat(Tz_hi_summer,365,1); 

Tz_lo_summer = [21*ones(8*12,1); 19*ones(10*12,1); 21*ones(6*12+1,1)];
Tz_lo_summer = repmat(Tz_lo_summer,365,1); 

% wall temperature setpoints
Tw_hi_summer = Tz_hi_summer + 5;
Tw_lo_summer = Tz_lo_summer - 5;


% Simulation w/ temperature deadband
for t = 1:control_dt/dt:length(SolWall) %1:length(solWall)
    if xx(2,n) > 24
        u(n) = -pHVAC;
    elseif xx(2,n) < 20
        u(n) = pHVAC;
    else
        u(n) = 0;
    end
    xx(1:2,n+1) = Ad*xx(1:2,n) + Bd*u(n) + Ed*W(:,t);
    xx(3:5,n) = W(:,t);
    n = n + 1;  % n is used to store x data, x and W have different time scale
end

% Plot results
figure; subplot(2,1,1); plot(xx(1:3, :)'); 
xticks([24*1*4 24*32*4 24*60*4 24*91*4 24*121*4 ...
    24*152*4 24*182*4 24*213*4 24*244*4 24*274*4 24*305*4 24*335*4 24*365*4])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
ylabel('Temperature (\circ C)')
legend({'Twall','Tzone','Tamb'}); title('Closed Loop Simulation with fixed HVAC power');
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec', 'Jan'})
subplot(2,1,2); plot(u); 
xticks([24*1*4 24*32*4 24*60*4 24*91*4 24*121*4 ...
    24*152*4 24*182*4 24*213*4 24*244*4 24*274*4 24*305*4 24*335*4 24*365*4])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec', 'Jan'})

% ON/OFF HEATING SEASON
% Simulation w/ temperature deadband
x0 = [22; 22];  % initial temperature values
x = zeros(5,length(SolWall)/(control_dt/dt)); % Twall, Tzone and Tamb
x(1:2,1) = x0;
heating_load_baseline = 0;
control_dt = 15*60;  % in seconds, SolWall is every dt second
n = 1; 

% Initialize states
xx = zeros(5,length(SolWall)/(control_dt/dt)); % Twall, Tzone and Tamb
xx(1:2,1) = x0; 


u = zeros(1,length(SolWall)/(control_dt/dt));
pHVAC = 4000;   %W
%n = 0;


for t = 6048:2:6048+2016 %1:length(solWall)
    if xx(2,n) > Tz_hi
        u(n) = -pHVAC;
    elseif xx(2,n) < Tz_lo
        u(n) = pHVAC;
        heating_load_baseline = heating_load_baseline + pHVAC;
    else
        u(n) = 0;
    end
    xx(1:2,n+1) = Ad*xx(1:2,n) + Bd*u(n) + Ed*W(:,t);
    xx(3:5,n) = W(:,t);
    n = n + 1;  % n is used to store x data, x and W have different time scale
end

% Plot results
figure; subplot(2,1,1); plot(xx(1:3, 1:2016/3)'); %(control_dt/dt)
xticklabels({'21', '22', '23', '24', '25', '26', '27'})
ylabel('Temperature (\circ C)')
legend({'Twall','Tzone','Tamb'}); title('Closed Loop Simulation with fixed HVAC power - Heating Season (January)');
subplot(2,1,2); plot(u(1:2016/(control_dt/dt))); 
xticklabels({'21', '22', '23', '24', '25', '26', '27'})
heating_load_baseline_final = heating_load_baseline * control_dt; % J

% ON/OFF COOLING SEASON
% Simulation w/ temperature deadband
%xx = xx(:,52417:54432);

xx = zeros(5,length(SolWall)/(control_dt/dt)); % Twall, Tzone and Tamb
xx(1:2,1) = x0; 


u = zeros(1,length(SolWall)/(control_dt/dt));
pHVAC = 3000;  % W

x0 = [22; 22];  % initial temperature values
xx = zeros(5,length(SolWall)/(control_dt/dt)); % Twall, Tzone and Tamb

cooling_load_baseline = 0;
heating_load_baseline = 0;
n = 52416/(control_dt/dt);

xx = zeros(5,length(SolWall)/(control_dt/dt)); % Twall, Tzone and Tamb
xx(1:2,n) = x0; 
for t = 52417+6048:3:54432+6048 %1:length(solWall) control_dt/dt
    if xx(2,n) > Tz_hi_summer
        u(n) = -pHVAC;
        cooling_load_baseline = cooling_load_baseline - pHVAC;
    elseif xx(2,n) < Tz_lo_summer
        u(n) = pHVAC;
        heating_load_baseline = heating_load_baseline + pHVAC;
    else
        u(n) = 0;
    end
    xx(1:2,n+1) = Ad*xx(1:2,n) + Bd*u(n) + Ed*W(:,t);
    xx(3:5,n) = W(:,t);
    n = n + 1;  % n is used to store x data, x and W have different time scale
end

cooling_load_baseline = (heating_load_baseline + abs(cooling_load_baseline))*control_dt*2;
%42416/3 --> 17466
%54432/3 --> 18144

% Plot results
figure; subplot(2,1,1); plot(xx(1:3,17466:18144)'); 
%plot(Tz_hi(52417+6048:52417+2016+6048), '--', 'color',[0.5,0.5,0.5])
%plot(Tz_lo(52417+6048:52417+2016+6048), '--', 'color',[0.5,0.5,0.5])
ylabel('Temperature (\circ C)')
xticklabels({'21', '22', '23', '24', '25', '26', '27'})
legend({'Twall','Tzone','Tamb'}); title('Closed Loop Simulation with fixed HVAC power Cooling Season (July)');
subplot(2,1,2); plot(u(17466:18144)); %19489:20160
xticklabels({'21', '22', '23', '24', '25', '26', '27'})

%% Import csv files
Weather = csvread('DryBulbTemperature(1).csv',2,0);
%Parameters = csvread('Parameters.csv',0,1);

%% Building parameters
%Additional Values
A_win = 2*3*2; %m^2
A_north = 2.7*8; %m^2
A_south = A_north - A_win;
A_east = 2.7*6;
A_west = 2.7*6;
A_wall = A_north + A_south + A_east + A_west;
A_roof = 8*6;
V_zone = 6*8*2.7;


R_wall_prime = 1.944 / A_wall;
R_roof_prime = 3.147 / A_roof;
R_inverse = 1/R_roof_prime + 1/R_wall_prime;
R = R_inverse^-1;


t_plasterboard = 0.012; %m
t_fiberglass = 0.066; %m
t_wood = 0.009; %m

rho_plasterboard = 950; %kg/m^3
rho_fiberglass = 12; %kg/m^3
rho_wood = 530; %kg/m^3

Cp_plasterboard = 840; %J/(kg*K)
Cp_fiberglass = 840; %J/(kg*K)
Cp_wood = 900; %J/(kg*K)

C_platerboard = rho_plasterboard*Cp_plasterboard*t_plasterboard*A_wall;
Cp_fiberglass = rho_fiberglass*Cp_fiberglass*t_fiberglass*A_wall;
C_wood = rho_wood*Cp_wood*t_wood*A_wall;

%Cw = C_platerboard + Cp_fiberglass + C_wood;  %J/K
R1 = R*0.5;
R2 = R1;

% Building parameters
%Wall_absorp = Parameters(1);
%Win_transm = Parameters(2);
%R1 = Parameters(3);
%R2 = Parameters(4);
%Rwin = Parameters(5);
%Cz = Parameters(6);
%Cw = Parameters(7);
%cop = Parameters(8);
%Qint_max = Parameters(9);
%Qsys = Parameters(10);

Wall_absorp = 0.6;
%Win_transm = win_transm(:);
%R1 = Handled above somehow
%R2 = Handled above somehow
%Rwin = 1/36;
%Cz = Handled above somehow
%Cw = Handled above somehow
%cop = 1;
%Qint_max = 200.* ones(1440,1);  %stated in problem --> Double Check
Qsys = 30;   % Double Check


% Building dimensions from ASHRAE Case 600
Asouth = 8*2.7-2*2*3;
Anorth = 8*2.7;
Aeast = 6*2.7;
Awest = 6*2.7;
Aroof = 6*8;
Awin = 2*2*3;

%% Build 3R2C model
A = [-(1/Cw)*(1/R1+1/R2), 1/(Cw*R2); ...
    1/(Cz*R2), -(1/Cz)*(1/R2+1/Rwin)];
B = [0; ...
    cop/Cz];
C = [0,1];
D = 0;
E = [1/(Cw*R1), 1/Cw, 0; ...
    1/(Cz*Rwin), 0, 1/Cz];

%% Discretize 
Ts = 60; % sampling time in seconds
Ad = expm(A*Ts);
Bd = inv(A)*(Ad-eye(2))*B;
Ed = inv(A)*(Ad-eye(2))*E;

%% Weather variables
%Qsol = Weather(:,1);
%Tamb = Weather(:,2);  % ambient temperatureTamb = Tamb_interp;  % ambient temperature
Tamb = Tamb_interp;  % ambient temperature
Qsol = (SolSouth*A_south+solNorth*A_north+SolEast*A_east+ ...
    SolWest*A_west+SolRoof*A_roof+SolWindow*A_win);

% plot(Tamb)
% plot(Qsol)

% Constant internal heat gains
Qint = 100*[ones(8*60*60/dt,1); zeros(10*60*60/dt,1); ones(6*60*60/dt,1)]; % out for work during the day from 8-18
% plot(Qint)
Qint = repmat(Qint,365,1);  % repmat([10;11],3,1); repeat matrix

%% temperature setpoints
% zone temperature 
Tz_hi = [24*ones(8*12,1); 26*ones(10*12,1); 24*ones(6*12+1,1)];
Tz_hi = repmat(Tz_hi,365,1); 

Tz_lo = [20*ones(8*12,1); 18*ones(10*12,1); 20*ones(6*12+1,1)];
Tz_lo = repmat(Tz_lo,365,1); 

% wall temperature setpoints
Tw_hi = Tz_hi + 5;
Tw_lo = Tz_lo - 5;

%% disturbance 

W = [Tamb_interp; SolWall'./100; Qint'];

%% Heating Baseline - Open loop simulation
x0 = [
    mean([Tw_lo(1),Tw_hi(1)]); 
    mean([Tz_lo(1),Tz_hi(1)])
    ];  % initial temperature values
x0 = [0;0];
x = zeros(2, length(W));
x(:,1) = x0;  % assign initial temperature values

u = zeros(1,length(W));
w = W; 

for t = 6048:6048+2016
    x(:,t+1) = Ad*x(:,t) +Bd*u(t) + Ed*w(:,t);
end

% PLOT
x =  x(:,6048:6048+2016);
figure
plot(x','LineWidth',1.25);
hold on
plot(Tamb(6048:6048+2016), 'LineWidth',1.25);
xticks([0 288*1, 288*2, 288*3, 288*4, 288*5, 288*6])
xticklabels({'21', '22', '23', '24', '25', '26', '27'});
ylabel('Temperature (\circ C)');ylim([0 40]);
legend({'Twall','Tzone','Tamb'}); 
title('Open Loop Simulation in Heating Season (January)');

%% Heating Baseline-  MPC 

Tp = 2*60*60/Ts;  % Prediction horizon = 2 hr
kb = 15*60/Ts;   % Control horizon = 15 min

x0 = [
    mean([Tw_lo(1),Tw_hi(1)]); 
    mean([Tz_lo(1),Tz_hi(1)])
    ];  % initial temperature values

% controllability matrices 
% A_bar 
A_bar = Ad;
for i = 2:Tp
    A_bar = [A_bar; Ad^i];
end

% B_bar
B_bar = zeros(2*Tp,Tp);
for i = 1:2:2*Tp
    B_bar(i:i+1,(i+1)/2) = Bd;
    xp = (i+1)/2 - 1;
    for j = 1:(i-1)/2
        B_bar(i:i+1,j) = Ad ^ xp * Bd;
        xp = xp - 1;
    end
end

% E_bar
E_bar = zeros(2*Tp,3*Tp);
for i = 1:2:2*Tp
    E_bar(i:i+1,(3*i-1)/2:(3*i+3)/2) = Ed;
    xp = (i+1)/2 - 1;
    for j = 1:3:(3*i-3)/2
        E_bar(i:i+1,j:j+2) = Ad ^ xp * Ed;
        xp = xp - 1;
    end
end

%% Heating Baseline - main algorithm - CVX
u_star = zeros(1,2016); % save the optimization results
x_star = zeros(2,2016);


x_star(:,1) = x0;

x_bar_hi = zeros(1, Tp*2);
x_bar_lo = zeros(1, Tp*2);

for t = 1:120:2016
    x_init = x_star(:,t);
    
    xx = zeros(2*(Tp+1), 1);
    
    % get the disturbance for this prediction window
    w_window = W(:,t+6048:t+6048+Tp-1);
    w_window = w_window(:);
    
    x_bar_hi(1:2:end) = Tw_hi(t+6048:t+6048+Tp-1);
    x_bar_hi(2:2:end) = Tz_hi(t+6048:t+6048+Tp-1);
    
    x_bar_lo(1:2:end) = Tw_lo(t+6048:t+6048+Tp-1);
    x_bar_lo(2:2:end) = Tz_lo(t+6048:t+6048+Tp-1);
    
    % call CVX
    cvx_begin
    cvx_solver sedumi
    % cvx_solver sdpt3
    variable u(Tp,1)
    minimize (norm(u,1)) % norm([1;2;3],1)
    
    subject to
        B_bar*u <= x_bar_hi' - A_bar*x_init - E_bar*w_window;
        -B_bar*u <= -x_bar_lo' + A_bar*x_init + E_bar*w_window;
        
        0 <= u <= 1e10*ones(Tp,1); % Heating
        %1e10*ones(Tp,1) <= u <= 0; % Heating
        
        u(1:kb) == repmat(u(1),[kb,1]);
        u(kb+1:2*kb) == repmat(u(kb+1),[kb,1]);
        u(2*kb+1:3*kb) == repmat(u(2*kb+1),[kb,1]);
        u(3*kb+1:4*kb) == repmat(u(3*kb+1),[kb,1]);
        u(4*kb+1:5*kb) == repmat(u(4*kb+1),[kb,1]);
        u(5*kb+1:6*kb) == repmat(u(5*kb+1),[kb,1]);
        u(6*kb+1:7*kb) == repmat(u(6*kb+1),[kb,1]);
        u(7*kb+1:8*kb) == repmat(u(7*kb+1),[kb,1]);
        
    cvx_end
    
    xx(1:2) = x_init;
    xx(3:end) = A_bar*x_init + B_bar*u + E_bar*w_window;
    
    x_star(:,t:t+Tp) = reshape(xx, [2,Tp+1]);
    u_star(t:t+Tp-1) = u;
        
end

%% Heating Baseline - plot figures 
figure
subplot(2,1,1);
plot(x_star(:,1:2016)', 'LineWidth',1.25)
hold on
plot(Tamb(6048:2016+6048),'LineWidth',1.25)
plot(Tz_hi(1:2016), '--', 'color',[0.5,0.5,0.5])
plot(Tz_lo(1:2016), '--', 'color',[0.5,0.5,0.5])
legend({'Twall','Tzone','Tamb'})
xticks([0 288*1, 288*2, 288*3, 288*4, 288*5, 288*6])
xticklabels({'21', '22', '23', '24', '25', '26', '27'});
ylabel('Temperature (\circ C)');ylim([0 40])
title('MPC Simulation - Heating Season (January)')

subplot(2,1,2)
plot(u_star, 'LineWidth',1.25)
xticks([0 288*1, 288*2, 288*3, 288*4, 288*5, 288*6])
xticklabels({'21', '22', '23', '24', '25', '26', '27'});
legend('pHVAC')
ylabel('HVAC Power (W)')
title('HVAC Power')

heating_mpc_load_sum = sum(u_star)*5*60;          % J
heating_mpc_load_max = max(u_star);

%% Cooling Baseline - Open loop simulation
x0 = [
    mean([Tw_lo(1),Tw_hi(1)]); 
    mean([Tz_lo(1),Tz_hi(1)])
    ];  % initial temperature values

x0 = [22;22];

x = zeros(2, length(w));
x = 20*ones(2, length(w));
x(:,1) = x0;  % assign initial temperature values

u = zeros(1,length(w));
w = W; 

for t = 52417+6048:52417+6048+2016
    x(:,t+1) = Ad*x(:,t) +Bd*u(t) + Ed*w(:,t);
end

x =  x(:,52417+6048:52417+6048+2016);
%x(:,1) = [22,22];
% PLOT
%close all
figure
plot(x','LineWidth',1.25);
hold on
plot(Tamb(52417+6048:52417+6048+2016), 'LineWidth',1.25);
xticks([0 288*1, 288*2, 288*3, 288*4, 288*5, 288*6])
xticklabels({'21', '22', '23', '24', '25', '26', '27'})
ylabel('Temperature (\circ C)');ylim([0 40]);
legend({'Twall','Tzone','Tamb'}); 
title('Open Loop Simulation in Cooling Season (July)');

%% Cooling Baseline-  MPC 

Tp = 2*60*60/Ts;  % Prediction horizon = 2 hr
kb = 15*60/Ts;   % Control horizon = 15 min

x0 = [
    mean([Tw_lo(1),Tw_hi(1)]); 
    mean([Tz_lo(1),Tz_hi(1)])
    ];  % initial temperature values

% controllability matrices 
% A_bar 
A_bar = Ad;
for i = 2:Tp
    A_bar = [A_bar; Ad^i];
end

% B_bar
B_bar = zeros(2*Tp,Tp);
for i = 1:2:2*Tp
    B_bar(i:i+1,(i+1)/2) = Bd;
    xp = (i+1)/2 - 1;
    for j = 1:(i-1)/2
        B_bar(i:i+1,j) = Ad ^ xp * Bd;
        xp = xp - 1;
    end
end

% E_bar
E_bar = zeros(2*Tp,3*Tp);
for i = 1:2:2*Tp
    E_bar(i:i+1,(3*i-1)/2:(3*i+3)/2) = Ed;
    xp = (i+1)/2 - 1;
    for j = 1:3:(3*i-3)/2
        E_bar(i:i+1,j:j+2) = Ad ^ xp * Ed;
        xp = xp - 1;
    end
end

%% Coling Baseline - main algorithm - CVX
u_star = zeros(1,2016); % save the optimization results
x_star = zeros(2,2016);


x_star(:,1) = x0;

x_bar_hi = zeros(1, Tp*2);
x_bar_lo = zeros(1, Tp*2);

for t = 1:24:2016
    x_init = x_star(:,t);
    
    xx = zeros(2*(Tp+1), 1);
    
    % get the disturbance for this prediction window
    w_window = W(:,t+52417+6048:t+52417+6048+Tp-1);
    w_window = w_window(:);
    
    x_bar_hi(1:2:end) = Tw_hi_summer(t+52417+6048:t+52417+6048+Tp-1);
    x_bar_hi(2:2:end) = Tz_hi_summer(t+52417+6048:t+52417+6048+Tp-1);
    
    x_bar_lo(1:2:end) = Tw_lo_summer(t+52417+6048:t+52417+6048+Tp-1);
    x_bar_lo(2:2:end) = Tz_lo_summer(t+52417+6048:t+52417+6048+Tp-1);
    
    % call CVX
    cvx_begin
    cvx_solver sedumi
    % cvx_solver sdpt3
    variable u(Tp,1)
    minimize (norm(u,1)) % norm([1;2;3],1)
    
    subject to
        B_bar*u <= x_bar_hi' - A_bar*x_init - E_bar*w_window;
        -B_bar*u <= -x_bar_lo' + A_bar*x_init + E_bar*w_window;
        
       0%1 <= u <= 1e10*ones(Tp,1); % Heating
        -1e6*ones(Tp,1) <= u <= 1e3*ones(Tp,1); % Cooling
        
        u(1:kb) == repmat(u(1),[kb,1]);
        u(kb+1:2*kb) == repmat(u(kb+1),[kb,1]);
        u(2*kb+1:3*kb) == repmat(u(2*kb+1),[kb,1]);
        u(3*kb+1:4*kb) == repmat(u(3*kb+1),[kb,1]);
        u(4*kb+1:5*kb) == repmat(u(4*kb+1),[kb,1]);
        u(5*kb+1:6*kb) == repmat(u(5*kb+1),[kb,1]);
        u(6*kb+1:7*kb) == repmat(u(6*kb+1),[kb,1]);
        u(7*kb+1:8*kb) == repmat(u(7*kb+1),[kb,1]);
        
    cvx_end
    
    xx(1:2) = x_init;
    xx(3:end) = A_bar*x_init + B_bar*u + E_bar*w_window;
    
    x_star(:,t:t+Tp) = reshape(xx, [2,Tp+1]);
    u_star(t:t+Tp-1) = u;
        
    
end

%% Cooling Baseline - plot figures 
figure
subplot(2,1,1);
plot(x_star(:,1:length(u_star))', 'LineWidth',1.25)
hold on
plot(Tamb(52417+6048:52417+2016+6048),'LineWidth',1.25)
plot(Tz_hi_summer(52417+6048:52417+2016+6048), '--', 'color',[0.5,0.5,0.5])
plot(Tz_lo_summer(52417+6048:52417+2016+6048), '--', 'color',[0.5,0.5,0.5])
legend({'Twall','Tzone','Tamb'})
xticks([0 288*1, 288*2, 288*3, 288*4, 288*5, 288*6])
xticklabels({'21', '22', '23', '24', '25', '26', '27'})
ylabel('Temperature (\circ C)');ylim([0 40])
title('MPC Simulation - Cooling Season (July)')

subplot(2,1,2)
plot(u_star, 'LineWidth',1.25)
xticks([0 288*1, 288*2, 288*3, 288*4, 288*5, 288*6])
xticklabels({'21', '22', '23', '24', '25', '26', '27'})
legend('pHVAC')
ylabel('HVAC Power (W)')
title('HVAC Power')

cooling_mpc_load_sum = sum(abs(u_star))*5*60;   % J
cooling_mpc_load_max = max(u_star);       % J
cooling_mpc_load_min = min(u_star);       % J

heating_percent_change = (heating_mpc_load_sum - heating_load_baseline_final)/heating_load_baseline_final;
cooling_percent_change = (cooling_mpc_load_sum - cooling_load_baseline)/cooling_load_baseline;

heating_maxpower_percent_change = (heating_mpc_load_max - 4000)/4000;
cooling_maxpower_percent_change = (-cooling_mpc_load_min - 3000)/3000;