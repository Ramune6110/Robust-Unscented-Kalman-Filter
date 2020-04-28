classdef  Adaptive_Robust_Unscented_Kalman_Filter < Autonomous_Mobile_Robot
    properties (Access = public)
        % ARUKF Weight Parameter
        alpha     = 0.001;       % Design parameter
        beta      = 2;           % Design parameter
        kappa     = 0;           % Design parameter
    end
    properties (Access = private)
        Omega = 0.15;   % The kernel width of Correntropy
        Lamda = 0.95;  % Lamda for Fading factor
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
        end
        function AMUKF = Adaptive_Robust_Unscented_kalman_filter(AMUKF, i)
            % Filtering Setup
            AMUKF.time                 = AMUKF.time + AMUKF.dt;
            AMUKF                      = Input(AMUKF);
            AMUKF                      = Observation(AMUKF);
            % Prediction Step
            AMUKF.Flag_ganerate_sigma  = 1;
            AMUKF                = GenerateSigmaPoints(AMUKF, AMUKF.Flag_ganerate_sigma);
            AMUKF                = PredictMotion(AMUKF);
            AMUKF.xPred                = (AMUKF.wm * AMUKF.sigma')';
            AMUKF.Flag_calculate_sigma = 1;
            AMUKF.PPred                = CalcSimgaPointsCovariance(AMUKF, AMUKF.Flag_calculate_sigma);
            % Filtering Step
            AMUKF.Flag_ganerate_sigma  = 0;
            AMUKF                = GenerateSigmaPoints(AMUKF, AMUKF.Flag_ganerate_sigma);
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
    end 
end