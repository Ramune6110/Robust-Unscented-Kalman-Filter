classdef  QS_Adaptive_Robust_Unscented_Kalman_Filter < Autonomous_Mobile_Robot & Common_Function
    properties (Access = public)
        % QS-ARUKF Weight Parameter
        alpha = 0.001;       % Design parameter
        beta  = 2;           % Design parameter
        kappa = 0;           % Design parameter
    end
    properties (Access = private)
        Omega = 0.15;   % The kernel width of Correntropy
        rho   = 0.95;  % Lamda for Fading factor
        ita   = 45;
    end
    properties (Access = private)
        zSigma               % sigma point mearsurment
        zb                   % Estimate Measurment value
        St                   % calculate sigma point
        Pxz                  % cross covariance matrix
        Pzz                  % measurment covariance matrix
        v                    % Residual 
        zeta                 % Residual
        phi                  % phi
        Rtilde               % Rtilda
        Lamda                % Lamda
        C_0                  % Fading factor for C
        Flag_ganerate_sigma  % Switch Prediction  1 or Filtering 0
        Flag_calculate_sigma % Switch Prediction  1 or Filtering 0 or Pzz 2
        Flag_outlier
    end
    methods (Access = public)
        % Initialize
        function QS_Adaptive_Robust_Unscented_Kalman_Filter = QS_Adaptive_Robust_Unscented_Kalman_Filter(xTrue, z, xEst, PEst, Q, R, Qsigma)
            QS_Adaptive_Robust_Unscented_Kalman_Filter@ Autonomous_Mobile_Robot(xTrue, z, xEst, PEst, Q, R, Qsigma);
        end
        function QS_ARUKF = QS_Adaptive_Robust_Unscented_kalman_filter(QS_ARUKF, i)
            % Filtering Setup
            QS_ARUKF.time                 = QS_ARUKF.time + QS_ARUKF.dt;
            QS_ARUKF                      = Input(QS_ARUKF);
            QS_ARUKF                      = Observation(QS_ARUKF);
            % Prediction Step
            QS_ARUKF.Flag_ganerate_sigma  = 1;
            QS_ARUKF                      = GenerateSigmaPoints(QS_ARUKF, QS_ARUKF.Flag_ganerate_sigma);
            QS_ARUKF                      = PredictMotion(QS_ARUKF);
            QS_ARUKF.xPred                = (QS_ARUKF.wm * QS_ARUKF.sigma')';
            QS_ARUKF.Flag_calculate_sigma = 1;
            QS_ARUKF.PPred                = CalcSimgaPointsCovariance(QS_ARUKF, QS_ARUKF.Flag_calculate_sigma);
            % Filtering Step
            QS_ARUKF.Flag_ganerate_sigma  = 0;
            QS_ARUKF                      = GenerateSigmaPoints(QS_ARUKF, QS_ARUKF.Flag_ganerate_sigma);
            QS_ARUKF.zSigma               = PredictObservation(QS_ARUKF);
            QS_ARUKF.zb                   = (QS_ARUKF.wm * QS_ARUKF.sigma')';
            QS_ARUKF.Pxz                  = CalcPxz(QS_ARUKF);
            QS_ARUKF.v                    = QS_ARUKF.z - QS_ARUKF.zb;
            QS_ARUKF.zeta                 = QS_ARUKF.R^(-1/2) * (QS_ARUKF.z - QS_ARUKF.zb);
            QS_ARUKF.phi                  = CalcPhi(QS_ARUKF);
            QS_ARUKF.Rtilde               = (QS_ARUKF.R')^(1/2) / QS_ARUKF.phi * (QS_ARUKF.R)^(1/2);
            if i == 1
                QS_ARUKF.C_0 = QS_ARUKF.v * QS_ARUKF.v';
            else
                QS_ARUKF.C_0 = (QS_ARUKF.rho * QS_ARUKF.C_0 + QS_ARUKF.v * QS_ARUKF.v') / (1 + QS_ARUKF.rho);
            end
            QS_ARUKF.Lamda                = CalcLamda(QS_ARUKF);
            QS_ARUKF.PPred                = QS_ARUKF.Lamda * QS_ARUKF.PPred;
            QS_ARUKF.Flag_calculate_sigma = 2;
            QS_ARUKF.Pzz                  = CalcSimgaPointsCovariance(QS_ARUKF, QS_ARUKF.Flag_calculate_sigma);
            QS_ARUKF.Pxz                  = CalcPxz(QS_ARUKF);
            QS_ARUKF.Flag_outlier         = 0;
            QS_ARUKF                      = Outlier(QS_ARUKF, QS_ARUKF.Flag_outlier);
            QS_ARUKF.K                    = QS_ARUKF.Pxz / QS_ARUKF.Pzz;
            QS_ARUKF.xEst                 = QS_ARUKF.xPred + QS_ARUKF.K * (QS_ARUKF.z - QS_ARUKF.zb);
            QS_ARUKF.PEst                 = QS_ARUKF.PPred - QS_ARUKF.K * QS_ARUKF.Pzz * QS_ARUKF.K';
            QS_ARUKF                      = Root_Mean_Square_Error(QS_ARUKF);
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
            elseif Flag == 0
                nSigma = length(this.zSigma(1, :));
                d      = this.zSigma - repmat(this.zb, 1, nSigma);
                P      = this.R;
                for i = 1 : nSigma
                    P = P + this.wc(i) * d(:, i) * d(:, i)';
                end
            elseif Flag == 2
                nSigma = length(this.zSigma(1, :));
                d      = this.zSigma - repmat(this.zb, 1, nSigma);
                for i = 1 : nSigma
                    P = this.wc(i) * d(:, i) * d(:, i)';
                end
                P = this.Lamda * P + this.Rtilde;
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
        function Phi = CalcPhi(this)
            First_content  = (1 / sqrt(2 * pi) * this.Omega^3) * (exp(-(this.zeta(1, 1)^2) / 2 * this.Omega^2));
            Second_content = (1 / sqrt(2 * pi) * this.Omega^3) * (exp(-(this.zeta(2, 1)^2) / 2 * this.Omega^2));
            Third_content  = (1 / sqrt(2 * pi) * this.Omega^3) * (exp(-(this.zeta(3, 1)^2) / 2 * this.Omega^2));
            Phi = diag([First_content Second_content Third_content]);
        end
        function Lamda_k = CalcLamda(this)
            H = eye(3);
            Lamda_0 = trace(this.C_0 - this.ita * this.R - H * this.Q * H') / trace(H * (this.PPred - this.Q) * H');
            
            if Lamda_0 > 1
                Lamda_k = Lamda_0;
            elseif Lamda_0 <= 1
                Lamda_k = 1;
            end
        end
    end 
end