% ワークスペースにx_opt、u_opt、K_optが展開されてるなら↓はコメントアウト
% load("SolverOut_rb.mat");
% load("RegulatorOut_rb.mat"); 

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

bn1 = 3/4;
an1 = 3/2*r_0/l_1; 
wn1 = 3/2*g/l_1; 
cn1 = 3*3*c_1/(m_1*l_1^2); 
bn2 = 3/4;
an2 = 3/2*r_0/l_2; 
wn2 = 3/2*g/l_2; 
cn2 = 3*c_2/(m_2*l_2^2); 

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
