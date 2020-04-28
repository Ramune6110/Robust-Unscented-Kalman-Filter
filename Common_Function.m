classdef Common_Function
    properties (Access = public)
        % Unscented kalman Filter factor
        sigma                % sigma point system
    end
    % Unscented Kalman Filter function
    methods (Access = public)
        function this = GenerateSigmaPoints(this, Flag)
            if Flag == 1
                this.sigma = this.xEst;
                Psqrt = sqrtm(this.PEst);
                m     = length(this.xEst);
                % Positive direction
                for ip = 1 : m
                    this.sigma = [this.sigma this.xEst + this.gamma * Psqrt(:, ip)];
                end
                % Negative direction
                for in = 1 : m
                    this.sigma = [this.sigma this.xEst - this.gamma * Psqrt(:, in)];
                end
            else
                this.sigma = this.xPred;
                Psqrt = sqrtm(this.PPred);
                m     = length(this.xPred);
                % Positive direction
                for ip = 1 : m
                    this.sigma = [this.sigma this.xPred + this.gamma * Psqrt(:, ip)];
                end
                % Negative direction
                for in = 1 : m
                    this.sigma = [this.sigma this.xPred - this.gamma * Psqrt(:, in)];
                end
            end
        end
        function this = PredictMotion(this)
            % Sigma Points predition with motion model
            for i = 1 : length(this.sigma(1, :))
                this.sigma(:, i) = Motion_model(this, i);
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
    end
end