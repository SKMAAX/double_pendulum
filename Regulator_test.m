%% レギュレータのみのテスト用コード

load("pm/SolverOut_pm.mat"); 
x0 = zeros(6, 1);
C = pm.Constants;
plant = pm.Plant(x0);

u_stab = 0;
x_traj = x_opt;
u_traj = u_opt;
x_eq = x1;
u_eq = u_stab;
x_traj(:,C.N+1) = x_eq;
u_traj(:,C.N+1) = u_eq;


% ---- ゲイン計算準備 --------------------------------------------------------
regulator = Regulator(@pm.Plant, C, x0);
nx = C.nx;          % 状態次元
N  = C.N+1;         % 離散地平長 ※mainの計算結果の都合上、k=Nではなくk=N+1でのxが平衡点となる
dt = C.dt;          % 積分ステップ ( >0 )

SR = C.SR;

% メモリ確保
Kt   = zeros(nx, N);      % 時変フィードバックゲイン
Pnow = SR;                  % 積分の初期値 P(t_f)=SR

% ---- 1. 有限時間最適制御ゲイン Kt(:,k) ---------------------------
% 末端時刻 (k = N+1) のゲイン
Kt(:,N) = regulator.findK(Pnow, x_traj(:,N), u_traj(:,N));

% 逆向きに DRE を積分しながらゲインを更新
for k = N-1:-1:1
    % Runge–Kutta 4 で 1 ステップ (負の dt で後退積分)
    Pnow = regulator.G_dt(Pnow, x_traj(:,k+1), u_traj(:,k+1), -dt);

    % 現在時刻 k のゲイン
    Kt(:,k) = regulator.findK(Pnow, x_traj(:,k), u_traj(:,k));
end

% ---- 2. 無限時間安定化ゲイン Kstab -------------------------------
%   「P を 0 から後退積分 → Kstab が収束するまで繰返し」
Pnow = zeros(nx);                           % 初期化 P(∞)=0
Kstab_prev = inf;                           % 収束判定用
maxIter = 5000;
tol     = 1e-12;

for iter = 1:maxIter
    % ゲイン計算
    Bt = regulator.plant.getB(x_eq);
    R = regulator.getR();
    
    Kstab = (R\Bt'*Pnow)';
    % 収束チェック
    if norm(Kstab - Kstab_prev, 'fro') < tol
        break
    end
    Kstab_prev = Kstab;

    % 後退積分 1 ステップ
    Pnow = regulator.G_dt(Pnow, x_eq, u_eq, -dt);
end
%---------------------------------

N_end = C.N2;

% 振り上げ後、シミュレーション終了時刻まで時変ゲインを延長
K_opt = [Kt Kstab*ones(1,N_end-N)];


%% 振り上げ安定化のシミュレーション

u = zeros(1,N_end);
du = zeros(1,N_end);
x = zeros(length(x0),N_end);
dx = zeros(length(x0),N_end);

x(:,1) = x0;
for i = 1:N_end
    dx(:,i) = x_opt(:,i) - x(:,i);
    du(i) = K_opt(:,i)'*dx(:,i);
    % if i <= C.N
        % du(i) = 0; % 振り上げ中のゲイン=0にするテスト
    % end

    u(i) = u_opt(i) + du(i); 

    x(:,i+1) = plant.Gp_dt(x(:,i), u(i), C.dt);
end
x(:,end) = [];

%% 結果の描画
plot_results(x, u, C.t2);