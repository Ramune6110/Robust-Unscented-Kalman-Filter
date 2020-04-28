classdef  Unscented_Kalman_Filter < Autonomous_Mobile_Robot
    properties (Access = public)
        % UKF Weight Parameter
        alpha     = 0.001;       % Design parameter
        beta      = 2;           % Design parameter
        kappa     = 0;           % Design parameter
    end
    properties (Access = private)
        zSigma               % sigma point mearsurment
        zb                   % Estimate Measurment value
        St                   % calculate sigma point
        Pxz                  % cross covariance matrix
        Flag_ganerate_sigma  % Switch Prediction  1 or Filtering 0
        Flag_calculate_sigma % Switch Prediction  1 or Filtering 0
        Flag_outlier  
    end
    methods (Access = public)
        % Initialize
        function Unscented_Kalman_Filter = Unscented_Kalman_Filter(xTrue, z, xEst, PEst, Q, R, Qsigma)
            Unscented_Kalman_Filter@ Autonomous_Mobile_Robot(xTrue, z, xEst, PEst, Q, R, Qsigma);
        end
        function UKF = Unscented_kalman_Filter(UKF)
            % Filtering Setup
            UKF.time                 = UKF.time + UKF.dt;
            UKF                      = Input(UKF);
            UKF                      = Observation(UKF);
            % Prediction Step
            UKF.Flag_ganerate_sigma  = 1;
            UKF                      = GenerateSigmaPoints(UKF, UKF.Flag_ganerate_sigma);
            UKF                      = PredictMotion(UKF);
            UKF.xPred                = (UKF.wm * UKF.sigma')';
            UKF.Flag_calculate_sigma = 1;
            UKF.PPred                = CalcSimgaPointsCovariance(UKF, UKF.Flag_calculate_sigma);
            % Filtering Step
            UKF.Flag_ganerate_sigma  = 0;
            UKF                      = GenerateSigmaPoints(UKF, UKF.Flag_ganerate_sigma);
            UKF.zSigma               = PredictObservation(UKF);
            UKF.zb                   = (UKF.wm * UKF.sigma')';
            UKF.Flag_calculate_sigma = 0;
            UKF.St                   = CalcSimgaPointsCovariance(UKF, UKF.Flag_calculate_sigma);
            UKF.Pxz                  = CalcPxz(UKF);
            UKF.K                    = UKF.Pxz / UKF.St;
            UKF.Flag_outlier         = 0;
            UKF                      = Outlier(UKF, UKF.Flag_outlier);
            UKF.xEst                 = UKF.xPred + UKF.K * (UKF.z - UKF.zb);
            UKF.PEst                 = UKF.PPred - UKF.K * UKF.St * UKF.K';
            UKF                      = Root_Mean_Square_Error(UKF);
        end 
    end
    methods (Access = protected)
        function P = CalcSimgaPointsCovariance(this, Flag)
            if Flag == 1
                nSigma = length(this.sigma(1, :));
                d      = this.sigma - repmat(this.xPred, 1, nSigma);
                P      = this.Q;
                for i = 1 : nSigma
                    P = P + this.wc(i) * d(:, i) * d(:, i)';
                end
            else
                nSigma = length(this.zSigma(1, :));
                d      = this.zSigma - repmat(this.zb, 1, nSigma);
                P      = this.R;
                for i = 1 : nSigma
                    P = P + this.wc(i) * d(:, i) * d(:, i)';
                end
            end
        end
        function Sigma = PredictObservation(this)
            % Sigma Points predition with observation model
            for i = 1 : length(this.sigma(1, :))
                Sigma(:, i) = Mearsure_model(this, i);
            end
        end    
        function P = CalcPxz(this)
            nSigma = length(this.sigma(1, :));
            dx     = this.sigma - repmat(this.xPred, 1, nSigma);
            dz     = this.zSigma - repmat(this.zb, 1, nSigma);
            P      = zeros(length(this.sigma(:, 1)));
            for i = 1 : nSigma
                P = P + this.wc(i) * dx(:, i) * dz(:, i)';
            end
        end
        function Z = Mearsure_model(this, i)
            x = this.sigma(:, i);
            C = [1 0 0;
                 0 1 0;
                 0 0 1];
            Z = C * x;
        end
    end 
end