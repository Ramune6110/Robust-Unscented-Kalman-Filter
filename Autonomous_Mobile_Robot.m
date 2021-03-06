classdef Autonomous_Mobile_Robot
    properties (Access = public)
        % Setup Time 
        time     = 0;   % Start Time
        endtime  = 60;  % End Time
        dt       = 0.1; % Sampling Time
        % Setup Input
        T        = 10;           % [sec]
        V        = 1.0;          % [m/s]
        yawrate  = 5 / 180 * pi; % [rad/s]
    end
    properties (Access = private)
        % Outlier value
        Outlier_z = [50; 50 ;0.01]; % Measurment Outlier
    end
    properties (Access = public)
        Q        % The Process Covariance matrix
        R        % The Measurment Covariance matrix
        Qsigma   % The Input Covariance matrix
        u        % Input
        xTrue    % State quantities
        z        % Measurment quantities
        xEst     % State estimation
        PEst     % The Covariance Matrix estimation
        xPred    % Prior State estimation
        PPred    % The Prior Covariance Matrix estimation
        K        % Kalman Gain    
        zEst     % Measurment estimation
        nsteps   % Iteration Steps
        Flag_sys % Swicth xTrue 1 or xEst 0
        Flag_mea % Switch xTrue 1 or xPred 0
        RMSE     % Root Mean Square Error
    end
    properties (Access = public)
        % Sigma Weight
        n                    % size of state vector
        lamda                % Design parameter
        wm                   % weight wm
        wc                   % weight wc
        gamma                % Design parameter
    end
    methods (Access = public)
        % Initialize
        function this = Autonomous_Mobile_Robot(xTrue, z, xEst, PEst, Q, R, Qsigma)
            switch nargin
                case 7
                    this.xTrue  = xTrue;
                    this.z      = z;
                    this.xEst   = xEst;
                    this.PEst   = PEst;
                    this.Q      = Q;
                    this.R      = R;
                    this.Qsigma = Qsigma;
                    this.nsteps = ceil((this.endtime - this.time) / this.dt);
                otherwise
                    disp('Error');
            end
            % Setup Weight
            this.n     = length(xEst);
            this.lamda = this.alpha^2 * (this.n + this.kappa) - this.n;
            this.gamma = sqrt(this.n + this.lamda);
            % Calculate the weight
            this.wm    = [this.lamda / (this.lamda + this.n)];
            this.wc    = [(this.lamda / (this.lamda + this.n)) + (1 - this.alpha^2 + this.beta)];
            for i = 1 : 2 * this.n
                this.wm = [this.wm 1 / (2 * (this.n + this.lamda))];
                this.wc = [this.wc 1 / (2 * (this.n + this.lamda))];
            end
        end
        function this = Input(this)
            this.u  = [this.V*(1-exp(-this.time/this.T)); 
                       this.yawrate * (1-exp(-this.time/this.T))]; % [m/s rad/s]
        end
        function this = Root_Mean_Square_Error(this)
            this.RMSE = [sqrt(mean(this.xTrue(1) - this.xEst(1))^2);
                         sqrt(mean(this.xTrue(2) - this.xEst(2))^2);
                         sqrt(mean(this.xTrue(3) - this.xEst(3))^2);];
        end
        function ShowErrorEllipse(this, type)
            % caluclate eig, eig_valus
            [eig_vec, eig_valus] = eig(this.PEst(1:2, 1:2));
            % eig comparizon
            if eig_valus(1, 1) >= eig_valus(2, 2)
                long_axis  = 3 * sqrt(eig_valus(1, 1));
                short_axis = 3 * sqrt(eig_valus(2, 2));
                angle      = atan2(eig_vec(1, 2), eig_vec(1, 1));
            else
                long_axis  = 3 * sqrt(eig_valus(2, 2));
                short_axis = 3 * sqrt(eig_valus(1, 1));
                angle      = atan2(eig_vec(2, 2), eig_vec(2, 1));
            end
            % make Ellipse
            t = 0:10:360;
            x = [long_axis * cosd(t); short_axis * sind(t)];
            if(angle < 0)
                angle = angle + 2*pi;
            end
            % Ellipse Rotation
            Rr = [cos(angle) sin(angle); -sin(angle) cos(angle)];
            x = Rr * x;
            switch type
                case 0
                    plot(x(1, :) + this.xEst(1, 1), x(2, :) + this.xEst(2, 1), '-.b', 'linewidth', 1.0);
                case 1
                    plot(x(1, :) + this.xEst(1, 1), x(2, :) + this.xEst(2, 1), '-.g', 'linewidth', 1.0);
                case 2
                    plot(x(1, :) + this.xEst(1, 1), x(2, :) + this.xEst(2, 1), '-.m', 'linewidth', 1.0);
                otherwise
                    disp('Nothing');
            end
            % Ellipse Vector(3sigma range)
            v = 3 * sqrt(eig_valus(1, 1)) * eig_vec(:, 1);
            quiver(this.xEst(1, 1), this.xEst(2, 1), v(1, 1), v(2, 1), 'red', 'linewidth', 1.5); hold on;
            v = 3 * sqrt(eig_valus(2, 2)) * eig_vec(:, 2);
            quiver(this.xEst(1, 1), this.xEst(2, 1), v(1, 1), v(2, 1), 'blue', 'linewidth', 1.5); hold on;
        end
    end
    methods (Access = protected)
        % System Model
        function x = f(this, Flag)
            if Flag == 1
                A = [1 0 0;
                     0 1 0;
                     0 0 1];
    
                B = [this.dt * cos(this.xTrue(3)) 0;
                     this.dt * sin(this.xTrue(3)) 0;
                     0 this.dt];

                x = A*this.xTrue + B * this.u;
            else
                A = [1 0 0;
                     0 1 0;
                     0 0 1];
    
                B = [this.dt * cos(this.xEst(3)) 0;
                     this.dt * sin(this.xEst(3)) 0;
                     0 this.dt];

                x = A*this.xEst + B * this.u;
            end
        end
        % Measurement Model
        function z = h(this, Flag)
            C = [1 0 0;
                 0 1 0;
                 0 0 1];
            if Flag == 1
                z = C * this.xTrue;
            else
                z = C * this.xPred;
            end
        end
        function this = Observation(this)
            this.Flag_sys = 1;                                                    % Swicth xTrue 1 or xEst 0
            this.xTrue    = f(this, this.Flag_sys);                               % Ground Truth
            this.u        = this.u + sqrtm(this.Qsigma) * randn(2, 1);            %add Process Noise�@�@
            this.Flag_mea = 1;                                                    % Switch xTrue 1 or xPred 0
            this.z        = h(this, this.Flag_mea) + sqrtm(this.R) * randn(3, 1); %Simulate Observation   
        end
        % Outlier
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