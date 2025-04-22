classdef aCS
    % Adaptive conditional sampling with Nataf transformation
    properties(Constant)
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % SedSav     Save the seeds or not (defalut:false)            logical           [1,1]
        % AcpOpt     Optimal acceptance rate (defalut:0.44)           double            [1,1]
        % NumIter    Number of iterations (defalut:5)                 double            [1,1]
        % Lambda     Adjustment factor (defalut:0.6)                  double            [1,1]
        % -----------------------------------------------------------------------------------
        SedSav = false;
        AcpOpt = 0.44;
        NumIte = 5;
        Lambda = 0.6;
    end

    properties
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % BayStc     A Bayesian structure                             /                 [1,1]
        % Ndim       Number of dimensions                             double            [1,1]
        % Nfun       Number of functions                              double            [1,1]
        % -----------------------------------------------------------------------------------
        BayStc
        Ndim
        Nfun
    end

    methods
        function obj = aCS(BayStc)
            % Initializatipon
            % ----------------------------------------------------------------------
            % SYNTAX:
            % obj = MCMC(Liklihood,Prior)
            % ----------------------------------------------------------------------
            % INPUTS:
            % Liklihood : Likelihood function
            % Prior     : Nataf based prior distribution
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % obj       : constructed class
            % ----------------------------------------------------------------------
            obj.BayStc = BayStc;
            obj.Nfun = BayStc.Nfun;
            obj.Ndim = BayStc.Ndim;
        end

        function [u,y,Ncal] = SamGen(obj,varargin)
            % Generate samples distributing as prior * likelihood
            % 1. With seeds, no "burn in"
            % 2. Without seeds, prior --> prior * likelihood
            % ----------------------------------------------------------------------
            % SYNTAX:
            % [u,y,Nfun] = SamGen(obj,Ns,Nc)
            % [u,y,Nfun] = SamGen(obj,Ns,x_sed,y_sed)
            % ----------------------------------------------------------------------
            % INPUTS:
            % obj   : class constructed
            % Ns    : length of Markov chains                                  [1,1]
            % Nc    : number of Markov chains                                  [1,1]
            % u_sed : initial samples                                      [Ndim,Nc]
            % y_sed : likelihood function of u_sed                         [Nfun,Nc]
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % u     : generated samples                                 [Ndim,Nc,Ns]
            % y     : output function values                            [Nfun,Nc,Ns]
            %         --likelihood function L
            %         --...
            %         --prior distribution p
            %         --...
            % Ncal  : number of likelihood calls                               [1,1]
            % ----------------------------------------------------------------------
            % REFERENCES:
            % ----------------------------------------------------------------------
            % [1].
            % ----------------------------------------------------------------------

            % initialization
            baystc = obj.BayStc;
            sedsav = obj.SedSav || nargin<=3;
            nDim = obj.Ndim;
            nFun = obj.Nfun;
            Nite = obj.NumIte;
            lambda = obj.Lambda;
            acpopt = obj.AcpOpt;
            if nargin<=3
                Nc = varargin{1};
                Ns = varargin{2};
                u_sed = randn(nDim,Nc);
                baystc = baystc.EvlY(u_sed);
                y_sed = [baystc.EvlLKF;baystc.EvlPDF;baystc.Y];
                Ncal = baystc.Ncal;
            else
                u_sed = varargin{1};
                y_sed = varargin{2};
                Ns = varargin{3};
                Nc = size(u_sed,2);
                Ncal = 0;
            end

            % allocate memory
            u = zeros(nDim,Nc,Ns);
            y = zeros(nFun,Nc,Ns);

            % generation of randn numbers
            du = randn(nDim,Nc,Ns);
            p = rand(1,Nc,Ns);

            % number of group the chains
            Ngrp = max(ceil(Nc/Nite),1);
            % index of start and end points of seeds
            ind_grp = min(repmat([1,Ngrp],Nite,1)+(0:Ngrp:(Nite-1)*Ngrp)',Nc);
            % standard deviation for each RV
            sigma0 = max(std(u_sed,0,2),eps);

            for iter = 1:Nite
                sigma = min(lambda*sigma0,1); % adaptive standard deviation
                rho = sqrt(1-sigma.^2); % correlation between successive samples
                Nacp = 0;
                dNcal = 0;
                ind_sam = ind_grp(iter,1):ind_grp(iter,2);
                for k=1:Ns
                    if k == 1 && sedsav
                        u(:,ind_sam,k) = u_sed(:,ind_sam);
                        y(:,ind_sam,k) = y_sed(:,ind_sam);
                    else
                        u(:,ind_sam,k) = rho.*u_sed(:,ind_sam) + sigma.*du(:,ind_sam,k);
                        baystc = baystc.EvlY(u(:,ind_sam,k));
                        y(:,ind_sam,k) = [baystc.EvlLKF;baystc.EvlPDF;baystc.Y];
                        ind_rej =  exp(y(1,ind_sam,k)-y_sed(1,ind_sam)) < p(1,ind_sam,k);
                        ind_sam_rej = ind_sam(ind_rej);
                        u(:,ind_sam_rej,k) = u_sed(:,ind_sam_rej);
                        y(:,ind_sam_rej,k) = y_sed(:,ind_sam_rej);
                        Nacp = Nacp + length(ind_sam)-length(ind_sam_rej);
                        dNcal = dNcal + baystc.Ncal;
                        u_sed(:,ind_sam) = u(:,ind_sam,k);
                        y_sed(:,ind_sam) = y(:,ind_sam,k);
                    end
                end
                % average acceptance rate
                acp = Nacp/dNcal;
                % new scaling parameter
                lambda = exp(log(lambda)+(acp-acpopt)/sqrt(iter));
                % update Nfun
                Ncal = Ncal+dNcal;
            end
        end

    end
end