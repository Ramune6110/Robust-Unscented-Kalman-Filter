classdef  Adaptive_Robust_Unscented_Kalman_Filter < Autonomous_Mobile_Robot
    properties (Access = private)
        % AMUKF Parameter
        alpha = 0.001; % Design parameter
        beta  = 2;     % Design parameter
        kappa = 0;     % Design parameter
        Omega = 0.15;   % The kernel width of Correntropy
        Lamda = 0.95;  % Lamda for Fading factor
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
        Alpha                % Alpha
        C_0                  % Fading factor for C
        Flag_ganerate_sigma  % Switch Prediction  1 or Filtering 0
        Flag_calculate_sigma % Switch Prediction  1 or Filtering 0 or Pzz 2
        Flag_outlier         
    end
    methods (Access = public)
        % Initialize
        function Adaptive_Robust_Unscented_Kalman_Filter = Adaptive_Robust_Unscented_Kalman_Filter(xTrue, z, xEst, PEst, Q, R, Qsigma)
            Adaptive_Robust_Unscented_Kalman_Filter@ Autonomous_Mobile_Robot(xTrue, z, xEst, PEst, Q, R, Qsigma);
            % Setup
            Adaptive_Robust_Unscented_Kalman_Filter.n     = length(xEst);
            Adaptive_Robust_Unscented_Kalman_Filter.lamda = Adaptive_Robust_Unscented_Kalman_Filter.alpha^2 * (Adaptive_Robust_Unscented_Kalman_Filter.n + ... 
                                                            Adaptive_Robust_Unscented_Kalman_Filter.kappa) - Adaptive_Robust_Unscented_Kalman_Filter.n;
            Adaptive_Robust_Unscented_Kalman_Filter.gamma = sqrt(Adaptive_Robust_Unscented_Kalman_Filter.n + Adaptive_Robust_Unscented_Kalman_Filter.lamda);
            % Calculate the weight
            Adaptive_Robust_Unscented_Kalman_Filter.wm    = [Adaptive_Robust_Unscented_Kalman_Filter.lamda / (Adaptive_Robust_Unscented_Kalman_Filter.lamda ...
                                                             + Adaptive_Robust_Unscented_Kalman_Filter.n)];
            Adaptive_Robust_Unscented_Kalman_Filter.wc    = [(Adaptive_Robust_Unscented_Kalman_Filter.lamda / (Adaptive_Robust_Unscented_Kalman_Filter.lamda ...
                                                             + Adaptive_Robust_Unscented_Kalman_Filter.n)) ...
                                                             + (1 - Adaptive_Robust_Unscented_Kalman_Filter.alpha^2 + Adaptive_Robust_Unscented_Kalman_Filter.beta)];
            for i = 1 : 2 * Adaptive_Robust_Unscented_Kalman_Filter.n
                Adaptive_Robust_Unscented_Kalman_Filter.wm = [Adaptive_Robust_Unscented_Kalman_Filter.wm ...
                                                1 / (2 * (Adaptive_Robust_Unscented_Kalman_Filter.n + Adaptive_Robust_Unscented_Kalman_Filter.lamda))];
                Adaptive_Robust_Unscented_Kalman_Filter.wc = [Adaptive_Robust_Unscented_Kalman_Filter.wc ...
                                                1 / (2 * (Adaptive_Robust_Unscented_Kalman_Filter.n + Adaptive_Robust_Unscented_Kalman_Filter.lamda))];
            end
        end
        function AMUKF = Adaptive_Robust_Unscented_kalman_filter(AMUKF, i)
            % Filtering Setup
            AMUKF.time                 = AMUKF.time + AMUKF.dt;
            AMUKF                      = Input(AMUKF);
            AMUKF                      = Observation(AMUKF);
            % Prediction Step
            AMUKF.Flag_ganerate_sigma  = 1;
            AMUKF.sigma                = GenerateSigmaPoints(AMUKF, AMUKF.Flag_ganerate_sigma);
            AMUKF.sigma                = PredictMotion(AMUKF);
            AMUKF.xPred                = (AMUKF.wm * AMUKF.sigma')';
            AMUKF.Flag_calculate_sigma = 1;
            AMUKF.PPred                = CalcSimgaPointsCovariance(AMUKF, AMUKF.Flag_calculate_sigma);
            % Filtering Step
            AMUKF.Flag_ganerate_sigma  = 0;
            AMUKF.sigma                = GenerateSigmaPoints(AMUKF, AMUKF.Flag_ganerate_sigma);
            AMUKF.zSigma               = PredictObservation(AMUKF);
            AMUKF.zb                   = (AMUKF.wm * AMUKF.sigma')';
            AMUKF.v                    = AMUKF.z - AMUKF.zb;
            AMUKF.zeta                 = AMUKF.R^(-1/2) * (AMUKF.zb - AMUKF.z);
            AMUKF.phi                  = CalcPhi(AMUKF);
            AMUKF.Rtilde               = (AMUKF.R')^(1/2) / AMUKF.phi * (AMUKF.R)^(1/2);
            AMUKF.Flag_calculate_sigma = 0;
            AMUKF.St                   = CalcSimgaPointsCovariance(AMUKF, AMUKF.Flag_calculate_sigma);
            if i == 1
                AMUKF.C_0 = AMUKF.v * AMUKF.v';
            else
                AMUKF.C_0 = (AMUKF.Lamda * AMUKF.C_0 + AMUKF.v * AMUKF.v') / (1 + AMUKF.Lamda);
            end
            AMUKF.Alpha                = CalcAlpha(AMUKF);
            AMUKF.PPred                = AMUKF.Alpha * AMUKF.PPred;
            AMUKF.Flag_calculate_sigma = 2;
            AMUKF.Pzz                  = CalcSimgaPointsCovariance(AMUKF, AMUKF.Flag_calculate_sigma);
            AMUKF.Pxz                  = CalcPxz(AMUKF);
            AMUKF.Flag_outlier         = 0;
            AMUKF                      = Outlier(AMUKF, AMUKF.Flag_outlier);
            AMUKF.K                    = AMUKF.Pxz / AMUKF.Pzz;
            AMUKF.xEst                 = AMUKF.xPred + AMUKF.K * (AMUKF.z - AMUKF.zb);
            AMUKF.PEst                 = AMUKF.PPred - AMUKF.K * AMUKF.Pzz * AMUKF.K';
            AMUKF                      = Root_Mean_Square_Error(AMUKF);
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
                P = this.Alpha * P + this.Rtilde;
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
            P = this.Alpha * P;
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
        function alpha_k = CalcAlpha(this)
            alpha_0 = trace(this.C_0 - this.R) / trace(this.St - this.R);
            
            if alpha_0 > 1
                alpha_k = alpha_0;
            elseif alpha_0 <= 1
                alpha_k = 1;
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