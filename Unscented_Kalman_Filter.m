classdef  Unscented_Kalman_Filter < Autonomous_Mobile_Robot
    properties (Access = private)
        % UKF Parameter
        alpha     = 0.001;       % Design parameter
        beta      = 2;           % Design parameter
        kappa     = 0;           % Design parameter
        Outlier_z = [50; 50 ;0.01]; % Measurment Outlier
        Outlier_x = [5; 5; 0.01]; % System model Outlier
    end
    properties (Access = private)
        n                    % size of state vector
        lamda                % Design parameter
        wm                   % weight wm
        wc                   % weight wc
        gamma                % Design parameter
        sigma                % sigma point system
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
            % Setup
            Unscented_Kalman_Filter.n     = length(xEst);
            Unscented_Kalman_Filter.lamda = Unscented_Kalman_Filter.alpha^2 * (Unscented_Kalman_Filter.n + ... 
                                            Unscented_Kalman_Filter.kappa) - Unscented_Kalman_Filter.n;
            Unscented_Kalman_Filter.gamma = sqrt(Unscented_Kalman_Filter.n + Unscented_Kalman_Filter.lamda);
            % Calculate the weight
            Unscented_Kalman_Filter.wm    = [Unscented_Kalman_Filter.lamda / (Unscented_Kalman_Filter.lamda + Unscented_Kalman_Filter.n)];
            Unscented_Kalman_Filter.wc    = [(Unscented_Kalman_Filter.lamda / (Unscented_Kalman_Filter.lamda + Unscented_Kalman_Filter.n)) + ...
                                             (1 - Unscented_Kalman_Filter.alpha^2 + Unscented_Kalman_Filter.beta)];
            for i = 1 : 2 * Unscented_Kalman_Filter.n
                Unscented_Kalman_Filter.wm = [Unscented_Kalman_Filter.wm ...
                                                1 / (2 * (Unscented_Kalman_Filter.n + Unscented_Kalman_Filter.lamda))];
                Unscented_Kalman_Filter.wc = [Unscented_Kalman_Filter.wc ...
                                                1 / (2 * (Unscented_Kalman_Filter.n + Unscented_Kalman_Filter.lamda))];
            end
        end
        function UKF = Unscented_kalman_Filter(UKF)
            % Filtering Setup
            UKF.time                 = UKF.time + UKF.dt;
            UKF                      = Input(UKF);
            UKF                      = Observation(UKF);
            % Prediction Step
            UKF.Flag_ganerate_sigma  = 1;
            UKF.sigma                = GenerateSigmaPoints(UKF, UKF.Flag_ganerate_sigma);
            UKF.sigma                = PredictMotion(UKF);
            UKF.xPred                = (UKF.wm * UKF.sigma')';
            UKF.Flag_calculate_sigma = 1;
            UKF.PPred                = CalcSimgaPointsCovariance(UKF, UKF.Flag_calculate_sigma);
            % Filtering Step
            UKF.Flag_ganerate_sigma  = 0;
            UKF.sigma                = GenerateSigmaPoints(UKF, UKF.Flag_ganerate_sigma);
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
        function Sigma = GenerateSigmaPoints(this, Flag)
            if Flag == 1
                Sigma = this.xEst;
                Psqrt = sqrtm(this.PEst);
                m     = length(this.xEst);
                % Positive direction
                for ip = 1 : m
                    Sigma = [Sigma this.xEst + this.gamma * Psqrt(:, ip)];
                end
                % Negative direction
                for in = 1 : m
                    Sigma = [Sigma this.xEst - this.gamma * Psqrt(:, in)];
                end
            else
                Sigma = this.xPred;
                Psqrt = sqrtm(this.PPred);
                m     = length(this.xPred);
                % Positive direction
                for ip = 1 : m
                    Sigma = [Sigma this.xPred + this.gamma * Psqrt(:, ip)];
                end
                % Negative direction
                for in = 1 : m
                    Sigma = [Sigma this.xPred - this.gamma * Psqrt(:, in)];
                end
            end
        end
        function Sigma = PredictMotion(this)
            % Sigma Points predition with motion model
            for i = 1 : length(this.sigma(1, :))
                Sigma(:, i) = Motion_model(this, i);
            end
        end
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
        function X = Motion_model(this, i)
            x = this.sigma(:, i);
            A = [1 0 0;
                 0 1 0;
                 0 0 1];
            B = [this.dt * cos(x(3)) 0;
                 this.dt * sin(x(3)) 0;
                 0 this.dt];

            X = A*x + B * this.u;
        end
        function Z = Mearsure_model(this, i)
            x = this.sigma(:, i);
            C = [1 0 0;
                 0 1 0;
                 0 0 1];
            Z = C * x;
        end
         function this = Outlier(this, flag)
            if flag == 0
                if this.time >= 30 && this.time < 31
                    this.z = this.z + this.Outlier_z;
                end
            elseif flag == 1
                if this.time >= 30 && this.time < 31
                    this.z = this.z + this.Outlier_z;
                end
                if this.time >= 10 && this.time < 11
                    this.xPred = this.xPred + this.Outlier_x;
                end
            end
        end
    end 
end