classdef NumSim
    % Numerical post-processing based on subset simulation samples
    % ----------------------------------------------------------------------
    % SYNTAX:
    %   obj = NumSim(in)
    %   [x_pos, y_pos, x_all, y_all] = obj.PosSim()
    %   [cv_z_hat, cv_z] = obj.EvlVar()
    %
    % INPUT:
    %   in : Output structure from SuS (Subset Simulation)
    %
    % OUTPUT (properties & methods):
    %   Wevd   : Weights for evidence estimation
    %   Wpos   : Weights for posterior distribution
    %   Wres   : Resampling weights between levels
    %   Z      : Final log-evidence estimate
    %   Zi     : Log-evidence estimates at each level
    %   PosSim : Resampled posterior samples
    %   EvlVar : Coefficient of variation (COV) estimates for Z
    %
    % Reference:
    %   Handley & Lemos (2019) - Physical Review D 100.2:023512
    % ----------------------------------------------------------------------

    properties
        X       % Cell array of samples at each level
        Y       % Cell array of function values at each level
        G       % Cell array of intermediate parameters
        H       % (Unused or reserved)
        C       % Normalization constant (log) at each level
        NumIte  % Number of iterations
        NumSam  % Number of samples per level
        SamSze  % Sample matrix size at each level [Nchains × Nsamp/chain]
        IntBay  % Intermediate Bayesian model (SuS_Bay_Nataf or similar)
        Wevd    % Weights for evidence estimation
        Wpos    % Weights for posterior distribution
        Wres    % Resampling weights
        Z       % Estimated log-evidence
        Zi      % Log-evidence per iteration
    end

    methods
        function obj = NumSim(in)
            % Constructor
            % Extract and normalize input from SuS simulation
            obj.X = in.x;
            obj.Y = in.y;
            obj.G = in.g;
            obj.NumIte = in.Nite;
            obj.NumSam = in.Nsam;
            obj.SamSze = in.Nsze;
            obj.IntBay = in.intBay;

            Nite = obj.NumIte;
            y = obj.Y;
            g = obj.G;
            intBay = obj.IntBay;

            % Normalize log-likelihoods across iterations
            w_res = cell(Nite,1);
            c = zeros(Nite,1);
            for ite = 1:Nite-1
                intBay.G = g{ite+1};
                w_res{ite} = ResRat(intBay, y{ite});
                c(ite+1) = c(ite) + logmean(w_res{ite},"all");
                y{ite+1}(1,:) = y{ite+1}(1,:) - c(ite+1);
            end
            obj.C = c;

            % Compute weights and evidence
            w_evd = WgtEvd(g, y);
            w_pos = WgtPos(g, y);
            z_i = zeros(Nite,1);
            for ite = 1:Nite
                z_i(ite) = c(ite) + logmean(w_evd{ite},"all");
            end
            z = logsum(z_i,"all");

            % Normalize posterior weights
            for ite = 1:Nite
                w_pos{ite} = w_pos{ite} - z;
            end

            obj.Wevd = w_evd;
            obj.Wpos = w_pos;
            obj.Wres = w_res;
            obj.Z = z;
            obj.Zi = z_i;
        end

        function [x_pos, y_pos, x_all, y_all] = PosSim(obj)
            % Resample posterior samples based on weight
            x = obj.X;
            y = obj.Y;
            w_pos = obj.Wpos;

            [x_all, y_all] = VabCmb(x, y);       % Flattened samples
            w_pos = VabCmb(w_pos);              % Flattened weights

            Npos = round(exp(2*logsum(w_pos,"all") - logsum(2*w_pos,"all")));
            Nall = size(y_all,2);
            w_pos_norm = exp(w_pos - max(w_pos,[],"all"));

            ind_pos = randsample(1:Nall, Npos, true, w_pos_norm);
            x_pos = x_all(:, ind_pos);
            y_pos = y_all(:, ind_pos);
        end

        function [cv_z_hat, cv_z] = EvlVar(obj)
            % Evaluate coefficient of variation (COV) for evidence estimation

            w_res = obj.Wres;
            w_evd = obj.Wevd;
            z_i = obj.Zi;
            z = obj.Z;
            Nite = obj.NumIte;

            % Estimate variance for conditional probs p (resampling weights)
            c_p_hat = zeros(Nite,1);
            c_p = zeros(Nite,1);
            gamma_p = zeros(Nite,1);
            for ite = 1:Nite-1
                [c_p_hat(ite), c_p(ite), gamma_p(ite)] = CorMod(w_res{ite});
            end

            % Estimate variance for h (evidence weights)
            c_h_hat = zeros(Nite,1);
            c_h = zeros(Nite,1);
            gamma_h = zeros(Nite,1);
            for ite = 1:Nite
                [c_h_hat(ite), c_h(ite), gamma_h(ite)] = CorMod(w_evd{ite});
            end

            % Cross correlation between h and p
            rho_hp = zeros(Nite-1,1);
            rho_fi = zeros(Nite-1,1);
            gamma_hp = zeros(Nite-1,1);
            for ite = 1:Nite-1
                [rho_fi(ite), gamma_hp(ite)] = CorMod_cross(w_evd{ite}, w_res{ite});
                rho_hp(ite) = rho_fi(ite)*(1+gamma_hp(ite))/sqrt(1+gamma_p(ite))/sqrt(1+gamma_h(ite));
            end

            % Covariance matrix for log-evidence (diagonal and off-diagonal terms)
            cov_z_hat = zeros(Nite, Nite);
            cov_z = zeros(Nite, Nite);
            for ite = 1:Nite
                z1z1 = 2 * z_i(ite);
                if ite == 1
                    sumc_hat = c_h_hat(ite)^2;
                    sumc = c_h(ite)^2;
                else
                    sumcp_hat = sum(c_p_hat(1:ite-1).^2, "all");
                    sumcp = sum(c_p(1:ite-1).^2, "all");
                    sumc_hat = c_h_hat(ite)^2 + sumcp_hat;
                    sumc = c_h(ite)^2 + sumcp;
                end
                cov_z_hat(ite, ite) = z1z1 + log(sumc_hat);
                cov_z(ite, ite) = z1z1 + log(sumc);

                if ite < Nite
                    sumcp_hat = sum(c_p_hat(1:ite-1).^2, "all");
                    sumcp = sum(c_p(1:ite-1).^2, "all");
                    cpch_hat = c_p_hat(ite)*c_h_hat(ite)*rho_hp(ite);
                    cpch = c_p(ite)*c_h(ite)*rho_fi(ite);
                    sumc_hat = sumcp_hat + cpch_hat;
                    sumc = sumcp + cpch;

                    z1z2 = z_i(ite) + z_i(ite+1:Nite);
                    cov_z_hat(ite, ite+1:Nite) = z1z2 + log(sumc_hat);
                    cov_z(ite, ite+1:Nite) = z1z2 + log(sumc);
                    cov_z_hat(ite+1:Nite, ite) = cov_z_hat(ite, ite+1:Nite);
                    cov_z(ite+1:Nite, ite) = cov_z(ite, ite+1:Nite);
                end
            end

            cov_z_hat(isnan(cov_z_hat)) = -inf;

            % Final coefficient of variation (COV)
            cv_z_hat = exp(logsum(cov_z_hat,"all")/2 - z);
            cv_z = exp(logsum(cov_z,"all")/2 - z);
        end
    end
end
