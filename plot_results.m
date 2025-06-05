function plot_results(x_opt, u_opt, Tend)
    %--- 時間軸生成 --------------------------------------------------------
    t = linspace(0, Tend, length(u_opt));

    figure;

    %% 1) 位置 (左軸) と 速度 (右軸)
    subplot(3,1,1);

    % ---- 左：位置 -------------------------------------------------------
    yyaxis left
    plot(t, x_opt(1,:), 'b-', 'LineWidth', 2); hold on;
    ylabel('Angle  \theta [rad]');

    % ---- 右：速度 -------------------------------------------------------
    yyaxis right
    plot(t, x_opt(2,:), 'r--', 'LineWidth', 1.5);
    ylabel('\omega_\theta [rad/s]');

    grid on;
    title('Position and Velocity');
    legend({'Angle (x1)','Angular Vel(x2)'}, ...
           'Location','best');

    %% 2) 角度 (左軸) と 角速度 (右軸)
    subplot(3,1,2);

    % ---- 左：角度 -------------------------------------------------------
    yyaxis left
    plot(t, x_opt(3,:), 'b-', 'LineWidth', 2); hold on;
    plot(t, x_opt(5,:), 'g-', 'LineWidth', 2);
    ylabel('Angle  \phi [rad]');

    % ---- 右：角速度 -----------------------------------------------------
    yyaxis right
    plot(t, x_opt(4,:), 'b--', 'LineWidth', 1.5);
    plot(t, x_opt(6,:), 'g--', 'LineWidth', 1.5);
    ylabel('\omega_\phi [rad/s]');

    grid on;
    title('Angle and Angular Velocity');
    legend({'Angle 1(x3)','Angle 2(x5)','Angular Vel 1(x4)','Angular Vel 2(x6)'}, ...
           'Location','best');

    %% 3) 入力トルク
    subplot(3,1,3);
    plot(t, u_opt, 'c', 'LineWidth', 2);
    xlabel('Time  [s]');
    ylabel('Torque  [Nm]');
    grid on;

    %% 全体タイトル
    sgtitle('Pendulum Simulation Results');
end