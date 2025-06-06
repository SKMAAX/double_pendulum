% rb.Constants.m
classdef Constants
    properties (Constant)
        % システムパラメータ
        pi = 3.141592653589793;
        g = 9.8;
        
        % システム次元
        nx = 6;  % 状態変数の数
        ni = 1;  % 入力の数
        ns = 3;  % サンプリング数
        
        % 時間パラメータ
        t0 = 0.0;
        t1 = 1.8;
        t2 = 10.0;
        dt = 0.01;
        N = 181;  % t1/dt + 1
        N2 = 1001;  % t2/dt + 1
        
        % 振子が剛体モデルの場合↓
        r_0 = 0.3; % アーム半径 [m]
        l_1 = 0.2; % 剛体振子1 長さ [m]
        l_2 = 0.1; % 剛体振子2 長さ [m]
        c_1 = 0.001; % 剛体振子1 減衰係数 [N·m·s/rad]
        c_2 = 0.001; % 剛体振子2 減衰係数 [N·m·s/rad]
        m_1 = 1; % 剛体振子1の重さ [kg]
        m_2 = 1; % 剛体振子2の重さ [kg] 
        
        % beta_1 =  rb.Constants.l_1^2 / 3; % 剛体振子1 イナーシャ÷重さ
        % beta_2 =  rb.Constants.l_2^2 / 3; % 剛体振子2 イナーシャ÷重さ

        % ノミナルプラントパラメータ
        bn1 = 3/4;
        an1 = 3/2*rb.Constants.r_0 / rb.Constants.l_1;
        wn1 = 3/2*rb.Constants.g   / rb.Constants.l_1;
        cn1 = 3*rb.Constants.c_1 / (rb.Constants.m_1 * rb.Constants.l_1^2);
        bn2 = 3/4;
        an2 = 3/2*rb.Constants.r_0 / rb.Constants.l_2;
        wn2 = 3/2*rb.Constants.g   / rb.Constants.l_2;
        cn2 = 3*rb.Constants.c_2 / (rb.Constants.m_2 * rb.Constants.l_2^2);

        % モデル化誤差を入れた、実プラントパラメータ
        bp1 = rb.Constants.bn1 * 0.95;
        ap1 = rb.Constants.an1 * 0.95;
        wp1 = rb.Constants.wn1 * 0.95;
        cp1 = rb.Constants.cn1 * 0.95;
        bp2 = rb.Constants.bn2 * 0.95;
        ap2 = rb.Constants.an2 * 0.95;
        wp2 = rb.Constants.wn2 * 0.95;
        cp2 = rb.Constants.cn2 * 0.95;

        % ソルバーパラメータ
        k_max = 10000;
        
        % コスト関数の重み
        s = [10000.0; 10000.0; 80000.0; 5000.0; 80000.0; 5000.0];
        q = [1.0; 1.0; 1.0; 5.0; 1.0; 5.0];
        r = 1.0;
        
        % レギュレータの重み
        % 
        % sr = [1; 1; 10000.0; 100.0; 10000.0; 100.0]*1;
        % Kstabを求めるときに得た行列Pの値をSRとすると、終端時刻でのK(t1)がKstabに一致する
 SR = [    4.98413091743749    11.920780504494   -481.65232351516  -55.9357622250252   338.075585524614   27.5410410002082;
    11.920780504494     55.01952440644  -2344.19099276645  -272.238069043762   1656.78198696679    134.96832002407;
   -481.65232351516  -2344.19099276644   990093.440574772   114983.680956854  -811503.427525593   -66107.083644134;
  -55.9357622250251  -272.238069043762   114983.680956854   13353.6003663652  -94243.6135066594  -7677.31870424791;
   338.075585524613   1656.78198696678  -811503.427525591  -94243.6135066592    684372.07928747   55750.3730392181;
   27.5410410002082    134.96832002407  -66107.0836441339   -7677.3187042479   55750.3730392181   4541.58477711923];


        qr = [1.0; 1.0; 10.0; 1.0; 10.0; 1.0]*10;
        rr = 1000.0;
    end
end