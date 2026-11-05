classdef SuS_Bay_Nataf
    % Subset simulation with Nataf transformation for Bayesian updating

    properties
        % LKF     : Likelihood function object (must implement EvlLKF and Pdis.U2X)
        % Nfun    : Number of outputs per sample (default: 4)
        % Ndim    : Dimension of parameter space
        % X       : Samples in standard normal space [Ndim × Nsam]
        % Y       : Function values [4 × Nsam]
        MethodProfile = 'author_latest'
        LKF
        Nfun = 4;
        Ndim
        Ncal
        X
        Y
    end

    properties
        % G       : Parameters for intermediate likelihood function
        %           [L_i, p_i, L_i1, p_i1, z_i1, z_tmp1, z_sum1]
        G = struct('parameters', [-inf, 0, -inf, 0, -inf, -inf, -inf], 'lambda', 0.6);
    end

    properties(Dependent)
        % FlgCvg  : Convergence flag
        FlgCvg
    end

    methods
        function spec = sampler_profile(~,family,name)
            if nargin<3||isempty(name),name='author_latest';end
            spec=method_sampler_profile(family,name);
        end
        function obj = SuS_Bay_Nataf(LKF, method_profile)
            if nargin < 2 || isempty(method_profile), method_profile='author_latest'; end
            % Constructor
            % INPUT:
            %   LKF : likelihood function object
            obj.LKF = LKF;
            obj.Ndim = LKF.Ndim;
            obj.MethodProfile=method_profile;method_sampler_profile('sus_bayes',method_profile);obj.G.p=0.1;
        end

        function obj = EvlY(obj, u)
            % Evaluate function values for input samples
            % INPUT:
            %   u : samples in standard normal space [Ndim × Nsam]
            % OUTPUT:
            %   obj.Y : [likelihood; prior] log-values
            x = U2X(obj, u);
            L = obj.LKF.EvlLKF(x);
            P = logGauss(u);
            obj.X = u;
            obj.Y = [L; P];
            obj.Ncal = size(u, 2);
        end

        function Li = EvlLKF(obj)
            % Evaluate intermediate likelihood function (log indicator)
            y = obj.Y;
            L = y(1,:);
            g = obj.G.parameters;
            L_i = g(3);
            if isinf(L_i)
                Li = zeros(size(L));
            else
                Li = log(double(L > L_i));
            end
        end

        function pi = EvlPDF(obj)
            % Evaluate intermediate prior (log prior)
            y = obj.Y;
            p = y(2,:);
            pi = p;
        end

        function [g, obj] = UpdObj(obj, y, p)
            % Update parameters for intermediate distributions
            % INPUT:
            %   y : function values [4 × Nsam]
            %   p : target conditional probability
            g = obj.G.parameters;
            L_i = g(3);
            p_i = g(4);
            z_i = g(5);
            z_tmp = g(6);
            z_sum = g(7);
            L = y(3,:);
            Nsam = size(L,2);
            L_sort = sort(L, "descend");
            Nsed = round(Nsam * p);
            if Nsed > 0
                L_i1 = (L_sort(Nsed) + L_sort(Nsed+1)) / 2;
            else
                L_i1 = L_sort(1);
            end
            p_i1 = p_i + log(sum(L > L_i1) / Nsam);

            z_i1 = p_i + logmean(min(logminus(L,L_i), logminus(L_i1,L_i)), "all");
            z_tmp1 = p_i + logmean(logminus(L,L_i), "all");
            z_sum1 = logplus(logplus(logminus(z_sum,z_tmp), z_i), z_tmp1);

            g = [L_i, p_i, L_i1, p_i1, z_i1, z_tmp1, z_sum1];
            obj.G.parameters = g;
            g = obj.G;
        end

        function [u_sed, y_sed] = SltSed(obj, u, y, Nsam)
            % Select seed samples for next iteration
            u = u(1:obj.Ndim,:);
            y = y(1:obj.Nfun,:);
            L = y(3,:);
            g = obj.G.parameters;
            L_i1 = g(3);
            Ns = round(Nsam ./ (1:Nsam));
            Nsed = round(Nsam ./ Ns);
            ind_sed = L > L_i1;
            x_sed_can = u(:, ind_sed);
            y_sed_can = y(:, ind_sed);
            Nsed_can = sum(ind_sed, "all");
            Nsed = Nsed(find(Nsed_can - Nsed >= 0, 1, "last"));
            if Nsed_can >= Nsed
                randind = randperm(Nsed_can, Nsed);
            elseif Nsed_can > 0
                randind = randi(Nsed_can, 1, Nsed);
            else
                randind = [];
            end
            u_sed = x_sed_can(:, randind);
            y_sed = y_sed_can(:, randind);
            if ~isempty(y_sed)
                y_sed = obj.UpdY(y_sed);
            end
        end

        function y = UpdY(obj, y)
            % Update function values by recomputing Li and pi
            y = y(3:end,:);
            obj.Y = y;
            Pi = obj.EvlPDF;
            Li = obj.EvlLKF;
            y = [Li; Pi; y];
        end

        function FlgCvg = get.FlgCvg(obj)
            % Compute convergence flag based on change in thresholds and evidence
            g = obj.G.parameters;
            L_i = g(1);
            L_i1 = g(3);
            z_tmp = g(6);
            z_sum = g(7);
            FlgCvg = exp(logminus(L_i1,L_i) - logplus(L_i1,L_i)) < 1e-1 && ...
                     exp(z_tmp - z_sum) < 1e-1;
        end

        function x = U2X(obj, u)
            % Transform from standard normal to physical space via Nataf
            x = obj.LKF.Pdis.U2X(u);
        end
    end

    methods
        function options=sampler_options(obj,kernel)
            spec=method_sampler_profile('sus_bayes',obj.MethodProfile);options=struct();
            if isfield(spec.options,kernel),options=spec.options.(kernel);end
        end
        function defaults=run_defaults(obj)
            defaults=struct('NumSam',1000,'ConPrb',0.1,'MaxIte',10000,'Verbose',true);
            if isprop(obj,'Nsam'),defaults.NumSam=obj.Nsam;end
            if isfield(obj.G,'p'),defaults.ConPrb=obj.G.p;end
        end
        function c=start_run(~,owner,sampler)
            owner.IntBay.G.p=owner.ConPrb;
            c=SeMCS.start_context(owner,sampler,max(1,owner.MaxIte-1),'overlay');
        end
        function [c,r]=sample_level(~,c,level)
            [c,r]=SeMCS.sample_context(c,level);
        end
        function c=advance_level(~,c,r,level)
            if level+1<=numel(c.owner.G)
                c.model.G=SeMCS.prescribed_state(c.model.G,c.owner.G{level+1},false,'overlay');extra=0;
            else
                [~,c.model]=c.model.UpdObj(r.y,c.model.G.p);extra=0;
            end
            c.Ncal(level)=c.Ncal(level)+extra;c.g{level+1,1}=c.model.G;
            c.sampler.BayStc=c.model;
        end
        function c=select_run_seeds(~,c,r,level)
            if c.targets(level+1)<=0,error('SeMCS:ScheduleExhausted','Sample budget exhausted before seed selection.');end
            target=c.targets(level+1);[c.xs,c.ys]=c.model.SltSed(r.x,c.y{level},target);
            c.nc=size(c.xs,2);c.steps=round(target/c.nc);
            if c.owner.Verbose
                fprintf('Level %d: N=%d, chains=%d, steps=%d, calls=%d\n',level,r.count,r.shape(1),r.shape(2),c.Ncal(level));
            end
        end
        function [yes,c]=run_complete(~,c,~,phase)
            yes=false;if ~strcmp(phase,'after_seeds'),return;end
            if c.nc<=0 || numel(c.x)==c.limit,yes=true;return;end
            c.complete=c.model.FlgCvg;yes=c.complete;
        end
        function out=finish_run(~,c)
            out=SeMCS.pack_output(c);
        end
    end
end

function spec = method_sampler_profile(family, name)
% Sampling settings for this method; thresholds, budgets and stopping are unchanged.
if nargin < 2 || isempty(name), name='author_latest'; end
family=char(family); name=char(name);
if ~any(strcmp(name,{'author_latest','paper'}))
    error('method_sampler_profile:Name','Choose author_latest or paper.');
end
if ~strcmp(family,'sus_bayes')
    error('method_sampler_profile:Family','Unsupported method family for this model.');
end
options=struct();
options.aCS=struct('Sigma0Mode','seed_std','CarryLambda',false);
options.ES=struct('ESMode','coordinate','Ntry',20);
options.aCS.AcceptArithmetic='log';options.aCS.CoordinateMode='identity';
options.aCS.GroupMode='available';options.aCS.AcceptanceRateMode='attempts';
options.ES.AcceptArithmetic='log';options.ES.CoordinateMode='identity';options.ES.SeedUpdateMode='each_step';
options.CW_aMH=struct('AcceptArithmetic','log','PriorMode','gaussian','CoordinateMode','identity');
if strcmp(name,'author_latest')
    options.aCS.AcceptArithmetic='exp';options.ES.AcceptArithmetic='exp';options.CW_aMH=struct('AcceptArithmetic','exp');
    options.aCS.GroupMode='fixed';
end
spec=struct('family',family,'name',name,'default_sampler','aCS','options',options);
end
