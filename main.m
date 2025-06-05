% main.m
clear; close all;

% 初期値の設定
x0 = zeros(6, 1);
x1 = [0; 0; pi; 0; pi; 0];


%% モデルの設定：質点型振子(PointMass)→pm、剛体振子(RigidBody)→rb
model = "rb";  % ← "pm" か "rb" を選ぶ

switch model
    case "pm"
        import pm.Constants          
        plant = pm.Plant( x0 );
        C = pm.Constants;
        u0 = ones(1, pm.Constants.N);
        solver = Solver(@pm.Plant, pm.Constants, x0, x1, u0);
        reglator = Regulator(@pm.Plant, pm.Constants, x0);
    case "rb"
        import rb.Constants
        plant = rb.Plant( x0 );
        C = rb.Constants;
        u0 = ones(1, rb.Constants.N);
        solver = Solver(@rb.Plant, rb.Constants, x0,  x1, u0);
        reglator = Regulator(@rb.Plant, rb.Constants, x0);
end

%% 振子振り上げを実現する制御入力uの設計
% ソルバーの実行 (数分かかる)
solver.Solve_Control_Problem();
% load("pm/SolverOut_pm.mat"); % すでに実行済みで、時短したい場合
% load("rb/SolverOut_rb.mat"); % すでに実行済みで、時短したい場合

% 最適解の取得
x_opt = solver.output_x();
u_opt = solver.output_u();
%% 振り上げ計画の結果描画
plot_results(x_opt, u_opt, C.t1);

%% 時変ゲイン、振り上げ後のレギュレータゲインの設計
u_stab = 0;
[Kt, Kstab] = reglator.Design(x_opt,u_opt,x1,u_stab);

Kt(:,end) = []; % 最後の時刻のゲインはKstabと同一なのでクリア

N_end = C.N2;

% 振り上げ後、シミュレーション終了時刻まで状態、入力、ゲインを時間延長
x_opt = [x_opt x1*ones(1,N_end-C.N)];
u_opt = [u_opt zeros(1,N_end-C.N)];
K_opt = [Kt Kstab*ones(1,N_end-C.N)];

%% 振り上げ安定化のシミュレーション
u = zeros(1,N_end);
du = zeros(1,N_end);
x = zeros(length(x0),N_end);
dx = zeros(length(x0),N_end);

x(:,1) = x0;
for i = 1:N_end
    dx(:,i) = x_opt(:,i) - x(:,i);
    du(i) = K_opt(:,i)'*dx(:,i);
    % du(i) = 0; % 時変ゲインなしでのテスト用
    u(i) = u_opt(i) + du(i); 

    x(:,i+1) = plant.Gp_dt(x(:,i), u(i), C.dt);
end
x(:,end) = [];

%% 振り上げ安定化結果の描画
plot_results(x, u, C.t2);
