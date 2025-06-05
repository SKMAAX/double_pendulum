classdef (Abstract) PlantBase < handle
    methods (Abstract)
        xdot = fn(obj,x,u);      % 名目モデル
        xdot = fp(obj,x,u);     % 実プラント
        xNext = Gn_dt(obj,x,u,dt);
        xNext = Gp_dt(obj,x,u,dt);
    end
end
