% -----------------------------------------------------------------------------
% File : Robust Unscented Kalman Filter
%
% Discription : State Estimation for Autonomous Mobile Robot 
%
% Environment : MATLAB_R2019a 
%
% Author : Ramune6110
%
% Copyright (c): 2020 Ramune6110
% -----------------------------------------------------------------------------
clear;
close all;
clc;
% Create subclasss instance
% x, z, xEst, PEst Q R Qsigma
UKF = Unscented_Kalman_Filter([0;0;0], [0;0;0], [0;0;0], eye(3), diag([0.1 0.1 toRadian(1)]).^2, ...
                              diag([1.5 1.5 toRadian(3)]).^2, diag([0.1 toRadian(20)]).^2);
AMUKF = Adaptive_Robust_Unscented_Kalman_Filter([0;0;0], [0;0;0], [0;0;0], eye(3), diag([0.1 0.1 toRadian(1)]).^2, ...
                                                 diag([1.5 1.5 toRadian(3)]).^2, diag([0.1 toRadian(20)]).^2);
QS_ARUKF = QS_Adaptive_Robust_Unscented_Kalman_Filter([0;0;0], [0;0;0], [0;0;0], eye(3), diag([0.1 0.1 toRadian(1)]).^2, ...
                                                           diag([1.5 1.5 toRadian(3)]).^2, diag([0.1 toRadian(20)]).^2);    

% Save data box
result.time          = [];
result.z             = [];
result.xTrue         = [];
result.xEst_UKF      = [];
result.RMSE_UKF      = [];
result.xEst_AMUKF    = [];
result.RMSE_AMUKF    = [];
result.xEst_QS_ARUKF = [];
result.RMSE_QS_ARUKF = [];

tic;
for i = 1 : UKF.nsteps
    % UKF
    UKF = Unscented_kalman_Filter(UKF);
    %AMUKF
    AMUKF = Adaptive_Robust_Unscented_kalman_filter(AMUKF, i);
    %QSAMUKF
    QS_ARUKF = QS_Adaptive_Robust_Unscented_kalman_filter(QS_ARUKF, i);
    
    %Animation
    Animation(UKF, AMUKF, QS_ARUKF, i)
    
    % Save data
    result.time          = [result.time; UKF.time];
    result.z             = [result.z; UKF.z'];
    result.xTrue         = [result.xTrue; UKF.xTrue'];
    result.xEst_UKF      = [result.xEst_UKF; UKF.xEst'];
    result.RMSE_UKF      = [result.RMSE_UKF; UKF.RMSE'];
    result.xEst_AMUKF    = [result.xEst_AMUKF; AMUKF.xEst'];
    result.RMSE_AMUKF    = [result.RMSE_AMUKF; AMUKF.RMSE'];
    result.xEst_QS_ARUKF = [result.xEst_QS_ARUKF; QS_ARUKF.xEst'];
    result.RMSE_QS_ARUKF = [result.RMSE_QS_ARUKF; QS_ARUKF.RMSE'];
end
toc;

% Drow a graph offline
DrowGraph(result);

function radian = toRadian(degree)
    radian = degree/180*pi;
end
