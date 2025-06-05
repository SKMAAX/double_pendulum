% Regulator.m
classdef Regulator
    properties (SetAccess = private)
        P
        Kt
        Kstab
        plant
        constants
    end
    
    methods
        function obj = Regulator(Plant,Constants,x0)
            obj.constants = Constants;
            obj.P = zeros(obj.constants.nx, obj.constants.nx);
            obj.Kt = zeros(obj.constants.ni, obj.constants.nx, obj.constants.N);
            obj.Kstab = zeros(obj.constants.ni, obj.constants.nx);
            obj.plant = Plant(x0);
        end
        function F = findK(obj, Pt, x, u)
                Bt = obj.plant.getB(x);
                R = obj.getR();
                F = (R\Bt'*Pt)';
        end
        function y = f(obj, Pt, x, u)
                At = obj.plant.getA(x, u);
                Bt = obj.plant.getB(x);
                Q = obj.getQ();
                R = obj.getR();
                % y = Pt*Bt*inv(R)*Bt'*Pt - Pt*At - At'*Pt - Q;
                y = Pt*Bt/R*Bt'*Pt - Pt*At - At'*Pt - Q;
        end
        function P = G_dt(obj, Pt, x, u, dt)
            w1 = 1.0/6.0;
            w2 = 1.0/3.0;
            w3 = 1.0/3.0;
            w4 = 1.0/6.0;
            
            dx_1 = obj.f(Pt, x, u) * dt;
            dx_2 = obj.f(Pt + dx_1 * 0.5, x, u) * dt;
            dx_3 = obj.f(Pt + dx_2 * 0.5, x, u) * dt;
            dx_4 = obj.f(Pt+ dx_3, x , u) * dt;
            
            P = Pt + (dx_1 * w1 + dx_2 * w2 + dx_3 * w3 + dx_4 * w4);
        end
        function Q = getQ(obj)
            Q = diag(obj.constants.qr);
        end
        
        function R = getR(obj)
            R = obj.constants.rr;
        end
        
        % x,uの軌道が与えられたとき、その軌道周りの時変最適レギュレータゲインを求める
        function [Kt, Kstab] = Design(obj, x_traj, u_traj, x_eq, u_eq)
            % ---- 事前準備 --------------------------------------------------------
            nx = obj.constants.nx;          % 状態次元
            ni = obj.constants.ni;          % 入力次元
            N  = obj.constants.N+1;           % 離散地平長 ※ x_traj(k=N)ではなく、x_traj(k=N+1)が平衡点x_eqとなる
            dt = obj.constants.dt;          % 積分ステップ ( >0 )
            x_traj(:,obj.constants.N+1) = x_eq;
            u_traj(:,obj.constants.N+1) = u_eq;
            % 終端重みの設定
            SR = obj.constants.SR;
        
            % メモリ確保
            Kt   = zeros(nx, N);        % 時変フィードバックゲイン
            Pnow = SR;                  % 積分の初期値 P(t_f)=SR
        
            % ---- 1. 有限時間最適制御ゲイン Kt(:,:,k) ---------------------------
            % 末端時刻 (k = N+1) のゲイン
            Kt(:,N) = obj.findK(Pnow, x_traj(:,N), u_traj(:,N));
        
            % 逆向きに DRE を積分しながらゲインを更新
            for k = N-1:-1:1
                % Runge–Kutta 4 で 1 ステップ (負の dt で後退積分)
                Pnow = obj.G_dt(Pnow, x_traj(:,k+1), u_traj(:,k+1), -dt);
        
                % 現在時刻 k のゲイン
                Kt(:,k) = obj.findK(Pnow, x_traj(:,k), u_traj(:,k));
            end
            obj.Kt = Kt;        % ← プロパティに保存
        
            % ---- 2. 無限時間安定化ゲイン Kstab -------------------------------
            %   「P を 0 から後退積分 → Kstab が収束するまで繰返し」
            Pnow = zeros(nx);                           % 初期化 P(∞)=0
            Kstab_prev = inf;                           % 収束判定用
            maxIter = 5000;
            tol     = 1e-12;
        
            for iter = 1:maxIter
                % ゲイン計算
                Kstab = obj.findK(Pnow, x_eq, u_eq);
        
                % 収束チェック
                if norm(Kstab - Kstab_prev, 'fro') < tol
                    break
                end
                Kstab_prev = Kstab;
        
                % 後退積分 1 ステップ
                Pnow = obj.G_dt(Pnow, x_eq, u_eq, -dt);
            end
            obj.Kstab = Kstab;         % プロパティに保存
        end
    end
end