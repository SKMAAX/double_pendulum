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
            y(4) = rb.Constants.bp1 * x(2)^2 * sin(x(3)) * cos(x(3)) - rb.Constants.ap1 * u(1) * cos(x(3)) - rb.Constants.wp1 * sin(x(3)) - rb.Constants.cp1 * x(4);
            y(5) = x(6);
            y(6) = rb.Constants.bp2 * x(2)^2 * sin(x(5)) * cos(x(5)) - rb.Constants.ap2 * u(1) * cos(x(5)) - rb.Constants.wp2 * sin(x(5)) - rb.Constants.cp2 * x(6);
        end

        % ノミナルプラントの状態方程式
        function y = fn(obj, x, u)
            y = zeros(6, 1);
            y(1) = x(2);
            y(2) = u(1);
            y(3) = x(4);
            y(4) = rb.Constants.bn1 * x(2)^2 * sin(x(3)) * cos(x(3)) - rb.Constants.an1 * u(1) * cos(x(3)) - rb.Constants.wn1 * sin(x(3)) - rb.Constants.cn1 * x(4);
            y(5) = x(6);
            y(6) = rb.Constants.bn2 * x(2)^2 * sin(x(5)) * cos(x(5)) - rb.Constants.an2 * u(1) * cos(x(5)) - rb.Constants.wn2 * sin(x(5)) - rb.Constants.cn2 * x(6);
        end

        function A = getA(obj, x, u)
            A = zeros(rb.Constants.nx, rb.Constants.nx);
            A(1,2) = 1;
            A(3,4) = 1;

            A(4,2) = 2 * rb.Constants.bn1 * x(2) * sin(x(3)) * cos(x(3));
            A(4,3) = rb.Constants.bn1 * x(2)^2 * cos(2*x(3)) + rb.Constants.an1 * u(1) * sin(x(3)) - rb.Constants.wn1 * cos(x(3));
            A(4,4) = -rb.Constants.cn1;
            
            A(5,6) = 1;
            A(6,2) = 2 * rb.Constants.bn2 * x(2) * sin(x(5)) * cos(x(5));
            A(6,5) = rb.Constants.bn2 * x(2)^2 * cos(2*x(5)) + rb.Constants.an2 * u(1) * sin(x(5)) - rb.Constants.wn2 * cos(x(5));
            A(6,6) = -rb.Constants.cn2;
        end

        function B = getB(obj, x)
            B = zeros(rb.Constants.nx, rb.Constants.ni);
            B(2,1) = 1;
            B(4,1) = -rb.Constants.an1 * cos(x(3));
            B(6,1) = -rb.Constants.an2 * cos(x(5));
        end
        
        function y = H_x(obj, x, x1, u, lam)
            Q = diag([rb.Constants.q(1) rb.Constants.q(2) rb.Constants.q(3) rb.Constants.q(4) rb.Constants.q(5) rb.Constants.q(6)]);
            temp = [0;
                lam(1)+2*rb.Constants.bn1*lam(4)*x(2)*sin(x(3))*cos(x(3))+2*rb.Constants.bn2*lam(6)*x(2)*sin(x(5))*cos(x(5));
                lam(4)*(rb.Constants.bn1*x(2)^2*cos(2*x(3)) + rb.Constants.an1*u(1)*sin(x(3)) - rb.Constants.wn1*cos(x(3)));
                lam(3) - rb.Constants.cn1*lam(4);
                lam(6)*(rb.Constants.bn2*x(2)^2*cos(2*x(5)) + rb.Constants.an2*u(1)*sin(x(5)) - rb.Constants.wn2*cos(x(5)));
                lam(5) - rb.Constants.cn2*lam(6);
                ]';
            y = ((x-x1)'*Q + temp)';
        end

        function y = H_u(obj, x, u, lam)
            y = rb.Constants.r * u(:) + lam(2) - lam(4) * rb.Constants.an1 * cos(x(3)) - lam(6) * rb.Constants.an2 * cos(x(5));
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
                obj.x = obj.Gn_dt(obj.x, u, rb.Constants.dt);
                t = t + rb.Constants.dt;
            end
        end

        function x_next = output(obj, x_t, u_t)
            x_next = obj.Gn_dt(x_t, u_t, rb.Constants.dt);
        end
    end
end