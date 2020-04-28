classdef  QS_Adaptive_Robust_Unscented_Kalman_Filter < Autonomous_Mobile_Robot
    properties (Access = private)
        % QS_ARUKF Parameter
        alpha = 0.001; % Design parameter
        beta  = 2;     % Design parameter
        kappa = 0;     % Design parameter
        Omega = 0.15;   % The kernel width of Correntropy
        rho   = 0.95;  % Lamda for Fading factor
        ita   = 45;
        Outlier_z = [50; 50; 0.01]; % Measurment Outlier
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
            % Setup
            QS_Adaptive_Robust_Unscented_Kalman_Filter.n     = length(xEst);
            QS_Adaptive_Robust_Unscented_Kalman_Filter.lamda = QS_Adaptive_Robust_Unscented_Kalman_Filter.alpha^2 * (QS_Adaptive_Robust_Unscented_Kalman_Filter.n + ... 
                                                            QS_Adaptive_Robust_Unscented_Kalman_Filter.kappa) - QS_Adaptive_Robust_Unscented_Kalman_Filter.n;
            QS_Adaptive_Robust_Unscented_Kalman_Filter.gamma = sqrt(QS_Adaptive_Robust_Unscented_Kalman_Filter.n + QS_Adaptive_Robust_Unscented_Kalman_Filter.lamda);
            % Calculate the weight
            QS_Adaptive_Robust_Unscented_Kalman_Filter.wm    = [QS_Adaptive_Robust_Unscented_Kalman_Filter.lamda / (QS_Adaptive_Robust_Unscented_Kalman_Filter.lamda ...
                                                             + QS_Adaptive_Robust_Unscented_Kalman_Filter.n)];
            QS_Adaptive_Robust_Unscented_Kalman_Filter.wc    = [(QS_Adaptive_Robust_Unscented_Kalman_Filter.lamda / (QS_Adaptive_Robust_Unscented_Kalman_Filter.lamda ...
                                                             + QS_Adaptive_Robust_Unscented_Kalman_Filter.n)) ...
                                                             + (1 - QS_Adaptive_Robust_Unscented_Kalman_Filter.alpha^2 + QS_Adaptive_Robust_Unscented_Kalman_Filter.beta)];
            for i = 1 : 2 * QS_Adaptive_Robust_Unscented_Kalman_Filter.n
                QS_Adaptive_Robust_Unscented_Kalman_Filter.wm = [QS_Adaptive_Robust_Unscented_Kalman_Filter.wm ...
                                                1 / (2 * (QS_Adaptive_Robust_Unscented_Kalman_Filter.n + QS_Adaptive_Robust_Unscented_Kalman_Filter.lamda))];
                QS_Adaptive_Robust_Unscented_Kalman_Filter.wc = [QS_Adaptive_Robust_Unscented_Kalman_Filter.wc ...
                                                1 / (2 * (QS_Adaptive_Robust_Unscented_Kalman_Filter.n + QS_Adaptive_Robust_Unscented_Kalman_Filter.lamda))];
            end
        end
        function QS_ARUKF = QS_Adaptive_Robust_Unscented_kalman_filter(QS_ARUKF, i)
            % Filtering Setup
            QS_ARUKF.time                 = QS_ARUKF.time + QS_ARUKF.dt;
            QS_ARUKF                      = Input(QS_ARUKF);
            QS_ARUKF                      = Observation(QS_ARUKF);
            % Prediction Step
            QS_ARUKF.Flag_ganerate_sigma  = 1;
            QS_ARUKF.sigma                = GenerateSigmaPoints(QS_ARUKF, QS_ARUKF.Flag_ganerate_sigma);
            QS_ARUKF.sigma                = PredictMotion(QS_ARUKF);
            QS_ARUKF.xPred                = (QS_ARUKF.wm * QS_ARUKF.sigma')';
            QS_ARUKF.Flag_calculate_sigma = 1;
            QS_ARUKF.PPred                = CalcSimgaPointsCovariance(QS_ARUKF, QS_ARUKF.Flag_calculate_sigma);
            % Filtering Step
            QS_ARUKF.Flag_ganerate_sigma  = 0;
            QS_ARUKF.sigma                = GenerateSigmaPoints(QS_ARUKF, QS_ARUKF.Flag_ganerate_sigma);
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
        function Sigma = GenerateSigmaPoints(this, Flag)
            if Flag == 1
                Sigma = this.xEst;
                [U, S, V] = svd(this.PEst);
                Psqrt = U * sqrt(S);
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
                [U, S, V] = svd(this.PPred);
                Psqrt = U * sqrt(S);
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
        function this = Outlier(this, flag)
            if flag == 0
                if this.time >= 30 && this.time <31
                    this.z = this.z + this.Outlier_z;
                end
            elseif flag == 1
                if this.time >= 30 && this.time <31
                    this.z = this.z + this.Outlier_z;
                end
                if this.time >= 10 && this.time < 11
                    this.xPred = this.xPred + this.Outlier_x;
                end
            end
        end
    end 
end