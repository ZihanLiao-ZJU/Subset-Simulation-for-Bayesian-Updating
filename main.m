classdef main
    % =======================================================================================
    %
    %       Sequential Multiple Importance Sampling Framework
    %       -------------------------------------------------
    %       Developed by:
    %           Zihan Liao, Zhejiang University
    %           Binbin Li, Zhejiang University
    %           Hua-Ping Wan, Zhejiang University
    %
    %       Initial Version: Jan 2023
    %       Last Modified  : Apr 2025
    %
    % =======================================================================================

    properties (Constant)
        % Maximum number of iterations
        MaxIte = 10000;
    end

    properties
        % Intermediate Bayesian model
        IntBay
        
        % Sampler object (must contain .BayStc, .SamGen methods)
        Sampler

        % Target conditional probability (e.g., for subset simulation)
        ConPrb = 0.1;

        % Number of samples per iteration (scalar or vector)
        NumSam = 1000;

        % Predefined intermediate distribution parameters
        G = [];

        % Storage for all samples
        X = [];
        Y = [];
    end

    methods
        function obj = main(Sampler)
            % Constructor for the main class
            % INPUT:
            %   Sampler : sampling object with BayStc, SamGen, etc.
            if nargin < 1
                error('Inputs must include the Sampler object with a valid BayStc field.');
            end
            obj.Sampler = Sampler;
            obj.IntBay = Sampler.BayStc;
        end

        function out = RunIte(obj)
            % Perform sequential sampling iterations
            %
            % OUTPUT:
            %   out : structure containing results of all iterations

            % Extract settings
            intBay  = obj.IntBay;
            sampler = obj.Sampler;
            p       = obj.ConPrb;
            maxi    = obj.MaxIte;
            flg_Nsam = numel(obj.NumSam) > 1;
            flg_IntD = ~isempty(obj.G);

            % Initialize storage
            Nsam = zeros(maxi,1);      % Sample count per iteration
            Nsze = zeros(maxi,2);      % [#chains, #samples per chain]
            Ncal = zeros(maxi,1);      % #function evaluations
            x = cell(maxi,1); y = cell(maxi,1); g = cell(maxi,1);

            % Sample size setup
            if flg_Nsam
                Nite_sam = length(obj.NumSam);
                Nsam(1:Nite_sam) = obj.NumSam;
            else
                Nsam(:) = obj.NumSam;
            end
            Nsze(1,:) = [Nsam(1), 1];

            % Intermediate distribution setup
            if flg_IntD
                Nite_int = length(obj.G);
                g(1:Nite_int) = obj.G;
            else
                Nite_int = 0;
                g{1} = intBay.G;
            end

            % Iterative sampling process
            for ite = 1:maxi
                % === Sample Generation ===
                if all(Nsze(ite,:) > 0) && ~isnan(Nsze(ite,2))
                    if ite == 1
                        [x{ite}, y{ite}, Ncal(ite)] = sampler.SamGen(Nsze(ite,1), Nsze(ite,2));
                    else
                        [x{ite}, y{ite}, Ncal(ite)] = sampler.SamGen(x_sed, y_sed, Nsze(ite,2));
                    end
                end
                Nsam(ite) = prod(Nsze(ite,:));

                % === Update Intermediate Distribution ===
                if ite + 1 <= Nite_int
                    intBay.G = g{ite+1};
                else
                    [g{ite+1}, intBay] = intBay.UpdObj(y{ite}, p);
                end
                sampler.BayStc = intBay;

                % === Seed Selection ===
                [x_sed, y_sed] = intBay.SltSed(x{ite}, y{ite}, Nsam(ite+1));
                Nsze(ite+1,1) = size(x_sed,2); % #chains
                Nsze(ite+1,2) = round(Nsam(ite+1) / Nsze(ite+1,1)); % samples/chain

                % === Display Summary ===
                DspRst(ite, Nsam(ite), Nsze(ite,:), Ncal(ite));

                % === Check for convergence or termination ===
                if Nsze(ite+1,1) <= 0 || ite + 1 == maxi
                    flg_cvg = false;
                    break
                elseif intBay.FlgCvg
                    flg_cvg = true;
                    break
                end
            end

            % === Truncate and Organize Outputs ===
            Nite = ite;
            [x, y, Nsam, Nsze, Ncal] = VabSrk(Nite, x, y, Nsam, Nsze, Ncal);
            g = VabSrk(Nite + 1, g);

            % === Output Struct ===
            out.x = x;
            out.y = y;
            out.g = g;
            out.Nsam = Nsam;
            out.Nsze = Nsze;
            out.Ncal = Ncal;
            out.Nite = Nite;
            out.flg_cvg = flg_cvg;
            out.intBay = intBay;
        end
    end
end
