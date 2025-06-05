classdef Solver < handle
    properties (SetAccess = private)
        x0
        x1
        x
        u
        lam %lambda
        g_J
        plant
        constants
    end

    methods
        function obj = Solver(Plant, Constants, x0, x1, u0)

            obj.plant = Plant(x0);
            obj.constants = Constants;
            if nargin == 5
                obj.x0 = x0;
                obj.x1 = x1;
                obj.u = u0;

                N  = numel(u0);           % 制御ステップ数
                nx = numel(x0);           % 状態次元
                ni = 1;                   % 入力数は1固定
                
                obj.x   = zeros(nx, N);
                obj.lam = zeros(nx, N);
                obj.g_J = zeros(ni, N);
            else
                error('Solver requires 5 input arguments: Plant, Constants, x0, x1, u0');
            end
        end


        % dx/dt = ∂H^T/∂λ, x(t_0) =x_0 の初期値問題を解く
        function IVP(obj)
            obj.x(:,1) = obj.x0;
            for i = 1:obj.constants.N-1
                obj.x(:,i+1) = obj.G_dt(obj.x(:,i), obj.u(:,i), obj.constants.dt);
            end
        end
        function x_next = G_dt(obj, x, u, dt)
            x_next = obj.plant.Gn_dt(x, u, dt);
        end


        % dλ/dt = ∂H^T/∂x, λ(t_1) = ∂Ψ^T/∂x の最終値問題を解く
        function FVP(obj)
            obj.lam(:,end) = -obj.constants.s.*(obj.x1 - obj.x(:,end));
            for i = obj.constants.N:-1:2
                obj.lam(:,i-1) = obj.G_dt2(obj.x(:,i), obj.u(:,i), obj.lam(:,i), -obj.constants.dt);
            end
        end
        function lam_next = G_dt2(obj, x, u, lam, dt)
            w = [1/6, 1/3, 1/3, 1/6];
            dlam_1 = obj.flam(x, u, lam) * dt;
            dlam_2 = obj.flam(x, u, lam + dlam_1 * 0.5) * dt;
            dlam_3 = obj.flam(x, u, lam + dlam_2 * 0.5) * dt;
            dlam_4 = obj.flam(x, u, lam + dlam_3) * dt;
            lam_next = lam + (dlam_1*w(1) + dlam_2*w(2) + dlam_3*w(3) + dlam_4*w(4));
        end
        function lam_next = flam(obj, x, u, lam)
            lam_next = -obj.plant.H_x(x, obj.x1, u,  lam);
        end

        % 評価関数の勾配方向を求める
        function g_J = grad_J(obj)
            g_J = zeros(obj.constants.ni, obj.constants.N);
            for i = 1:obj.constants.N
                g_J(:,i) = obj.plant.H_u(obj.x(:,i), obj.u(:,i), obj.lam(:,i));
            end
        end
        
        % 一次元(直線)探索アルゴリズムによる移動距離の算出
        function alpha = search_a(obj, u, dk)
            e = 1e-6;
            a1 = 0.0;
            h = 0.01;
            a2 = a1 + h;

            while obj.J(obj.u + dk*a1) < obj.J(obj.u + dk*a2)
                h = h / 10;
                a2 = a1 + h;
            end

            a3 = a2 + 2*h;

            i=2;
            while obj.J(obj.u + dk*a2) > obj.J(obj.u + dk*a3)
                a1 = a2;
                a2 = a3;
                a3 = a2 + 2*i*h;
                i = i+1;
            end

            A = [a1^2, a1, 1; a2^2, a2, 1; a3^2, a3, 1];
            b = [obj.J(obj.u + dk*a1); obj.J(obj.u + dk*a2); obj.J(obj.u + dk*a3)];

            x = A\b;
            alpha_opt = -x(2)/(2*x(1));

            if alpha_opt < 0 || abs(a2 - a1) < e || abs(a3 - a2) < e
                alpha_opt = a2;
            end
            alpha = alpha_opt;
        end

        function J_val = J(obj, u)
            
            
            % obj.u = u;
            % obj.IVP(); 
            obj.x(:,1) = obj.x0;
            for i = 1:obj.constants.N-1
                obj.x(:,i+1) = obj.G_dt(obj.x(:,i), u(:,i), obj.constants.dt);
            end

            terminal_cost = 0.5 * (obj.x(:,end)-obj.x1)'*diag(obj.constants.s)*(obj.x(:,end)-obj.x1);
            % integral_cost = 0;
            % for i = 1:obj.constants.N
            %     integral_cost = integral_cost + 0.5*(obj.x(:,i)-obj.x1)'*diag(obj.constants.q)*(obj.x(:,i)-obj.x1) + 0.5*obj.constants.r*u(:,i)^2;
            % end
            L = @(xi,ui) 0.5*((xi-obj.x1).'*diag(obj.constants.q)*(xi-obj.x1) + ...
                  obj.constants.r*ui(1)^2);

            integral_cost = 0.5*( L(obj.x(:,1),  u(:,1)) + ...
                                  L(obj.x(:,end),u(:,end)) );
            
            for i = 2:obj.constants.N-1
                integral_cost = integral_cost + L(obj.x(:,i), u(:,i));
            end
            integral_cost = integral_cost * obj.constants.dt;
            J_val = terminal_cost + integral_cost;
            
            % disp(['x(1): ',num2str(obj.x(:,1)')]);
            % disp(['x(end): ',num2str(obj.x(:,end)')]);
            % disp(['phi: ',num2str(terminal_cost)]);
            % disp(['Integral_L: ',num2str(integral_cost)]);
        end

        function Solve_Control_Problem(obj)

            for k = 1:obj.constants.k_max
                obj.IVP();
                obj.FVP();
                dk = -obj.grad_J();
                alpha = obj.search_a(obj.u, dk);
                obj.u = obj.u + alpha*dk;
                
                if mod(k,100) == 0
                    % disp(['d_k: ', num2str(dk)]);
                    % disp(['alpha_k: ', num2str(alpha)]);
                    disp(['k:',num2str(k),'(kmax=',num2str(obj.constants.k_max),...
                        '), a*d_k:',num2str(norm(alpha*dk)),', J(u): ', num2str(obj.J(obj.u))]);
                end
                if norm(alpha*dk) < 1e-6
                    disp(['Converged at iteration: ', num2str(k)]);
                    break;
                end
            end
        end

        function u_opt = output_u(obj)
            u_opt = obj.u;
        end

        function x_opt = output_x(obj)
            x_opt = obj.x;
        end
    end
end