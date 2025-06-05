% Plant.m
classdef Plant < PlantBase 
    properties
        x
    end

    methods
        function obj = Plant(x0)
            obj.x = x0;
        end
        
        % 実プラントの状態方程式
        function y = fp(obj, x, u)
            y = zeros(6, 1);
            y(1) = x(2);
            y(2) = u(1);
            y(3) = x(4);
            y(4) = x(2)^2 * sin(x(3)) * cos(x(3)) - pm.Constants.ap1 * u(1) * cos(x(3)) - pm.Constants.wp1 * sin(x(3)) - pm.Constants.cp1 * x(4);
            y(5) = x(6);
            y(6) = x(2)^2 * sin(x(5)) * cos(x(5)) - pm.Constants.ap2 * u(1) * cos(x(5)) - pm.Constants.wp2 * sin(x(5)) - pm.Constants.cp2 * x(6);
        end

        % ノミナルプラントの状態方程式
        function y = fn(obj, x, u)
            y = zeros(6, 1);
            y(1) = x(2);
            y(2) = u(1);
            y(3) = x(4);
            y(4) = x(2)^2 * sin(x(3)) * cos(x(3)) - pm.Constants.an1 * u(1) * cos(x(3)) - pm.Constants.wn1 * sin(x(3)) - pm.Constants.cn1 * x(4);
            y(5) = x(6);
            y(6) = x(2)^2 * sin(x(5)) * cos(x(5)) - pm.Constants.an2 * u(1) * cos(x(5)) - pm.Constants.wn2 * sin(x(5)) - pm.Constants.cn2 * x(6);
        end

        function A = getA(obj, x, u)
            A = zeros(pm.Constants.nx, pm.Constants.nx);
            A(1,2) = 1;
            A(3,4) = 1;

            A(4,2) = 2 * x(2) * sin(x(3)) * cos(x(3));
            A(4,3) = x(2)^2 * cos(2*x(3)) + pm.Constants.an1 * u(1) * sin(x(3)) - pm.Constants.wn1 * cos(x(3));
            A(4,4) = -pm.Constants.cn1;
            
            A(5,6) = 1;
            A(6,2) = 2 * x(2) * sin(x(5)) * cos(x(5));
            A(6,5) = x(2)^2 * cos(2*x(5)) + pm.Constants.an2 * u(1) * sin(x(5)) - pm.Constants.wn2 * cos(x(5));
            A(6,6) = -pm.Constants.cn2;
        end
        
        function B = getB(obj, x)
            B = zeros(pm.Constants.nx, pm.Constants.ni);
            B(2,1) = 1;
            B(4,1) = -pm.Constants.an1 * cos(x(3));
            B(6,1) = -pm.Constants.an2 * cos(x(5));
        end
        
        function y = H_x(obj, x, x1, u, lam)
            y = zeros(pm.Constants.nx, 1);
            Q = diag([pm.Constants.q(1) pm.Constants.q(2) pm.Constants.q(3) pm.Constants.q(4) pm.Constants.q(5) pm.Constants.q(6)]);
            temp = [0;
                lam(1)+2*lam(4)*x(2)*sin(x(3))*cos(x(3))+2*lam(6)*x(2)*sin(x(5))*cos(x(5));
                lam(4)*(x(2)^2*cos(2*x(3)) + pm.Constants.an1*u(1)*sin(x(3)) - pm.Constants.wn1*cos(x(3)));
                lam(3) - pm.Constants.cn1*lam(4);
                lam(6)*(x(2)^2*cos(2*x(5)) + pm.Constants.an2*u(1)*sin(x(5)) - pm.Constants.wn2*cos(x(5)));
                lam(5) - pm.Constants.cn2*lam(6);
                ]';
            y = ((x-x1)'*Q + temp)';
        end

        function y = H_u(obj, x, u, lam)
            y = pm.Constants.r * u(:) + lam(2) - lam(4) * pm.Constants.an1 * cos(x(3)) - lam(6) * pm.Constants.an2 * cos(x(5));
        end

        % ノミナルプラントの時間発展
        function x_next = Gn_dt(obj, x_t, u_t, dt_t)
            w1 = 1.0/6.0;
            w2 = 1.0/3.0;
            w3 = 1.0/3.0;
            w4 = 1.0/6.0;

            dx_1 = obj.fn(x_t, u_t) * dt_t;
            dx_2 = obj.fn(x_t + dx_1 * 0.5, u_t) * dt_t;
            dx_3 = obj.fn(x_t + dx_2 * 0.5, u_t) * dt_t;
            dx_4 = obj.fn(x_t + dx_3, u_t) * dt_t;

            x_next = x_t + (dx_1 * w1 + dx_2 * w2 + dx_3 * w3 + dx_4 * w4);
        end

        % 実プラントの時間発展
        function x_next = Gp_dt(obj, x_t, u_t, dt_t)
            w1 = 1.0/6.0;
            w2 = 1.0/3.0;
            w3 = 1.0/3.0;
            w4 = 1.0/6.0;

            dx_1 = obj.fp(x_t, u_t) * dt_t;
            dx_2 = obj.fp(x_t + dx_1 * 0.5, u_t) * dt_t;
            dx_3 = obj.fp(x_t + dx_2 * 0.5, u_t) * dt_t;
            dx_4 = obj.fp(x_t + dx_3, u_t) * dt_t;

            x_next = x_t + (dx_1 * w1 + dx_2 * w2 + dx_3 * w3 + dx_4 * w4);
        end

        function Move(obj, t_1)
            t = 0.0;
            u = 1.0;

            while t <= t_1
                obj.x = obj.Gn_dt(obj.x, u, pm.Constants.dt);
                t = t + pm.Constants.dt;
            end
        end

        function x_next = output(obj, x_t, u_t)
            x_next = obj.Gn_dt(x_t, u_t, pm.Constants.dt);
        end
    end
end