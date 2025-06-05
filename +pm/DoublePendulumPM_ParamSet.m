% ワークスペースにx_opt、u_opt、K_optが展開されてるなら↓はコメントアウト
load("pm/SolverOut_pm.mat");
load("pm/RegulatorOut_pm.mat"); 

g = 9.80665;
dt = 0.01;
Tend = 10;

r_0 = 0.3;
l_1 = 0.2;
l_2 = 0.1;
c_1 = 0.001;
c_2 = 0.001;
m_1 = 1; % [kg]
m_2 = 1; % [kg] 


an1 = r_0/l_1; 
wn1 = g/l_1; 
cn1 = c_1/(m_1*l_1^2); 
an2 = r_0/l_2; 
wn2 = g/l_2; 
cn2 = c_2/(m_2*l_2^2); 
 

% x,u,Kの時間構造体を作成
x_opt_ts.time = linspace(0, Tend, length(u_opt));
x_opt_ts.signals.values = x_opt.';
x_opt_ts.dimensions = 6;

u_opt_ts.time = x_opt_ts.time;
u_opt_ts.signals.values = u_opt.';
u_opt_ts.dimensions = 1;

K_ts.time = x_opt_ts.time;
K_ts.signals.values = K_opt.';
K_ts.dimensions = 6;
