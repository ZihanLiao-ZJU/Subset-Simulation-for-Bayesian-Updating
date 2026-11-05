classdef SuS_Bay_Nataf
    % Subset simulation with Nataf transformation for Bayesian updating

    properties
        % LKF     : Likelihood function object (must implement EvlLKF and Pdis.U2X)
        % Nfun    : Number of outputs per sample (default: 4)
        % Ndim    : Dimension of parameter space
        % X       : Samples in standard normal space [Ndim × Nsam]
        % Y       : Function values [4 × Nsam]
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
        function spec = sampler_profile(obj,family,name)
            if nargin<3||isempty(name),name='author_latest';end
            if strcmp(name,'closed_ties_unit_fallback_v1')&&isfield(obj.G,'method_profile')&& ...
                    ~strcmp(obj.G.method_profile,'closed_ties_unit_fallback_v1')
                error('SuSClosedFallback:FreshModelRequired','Construct a fresh model with the closed-boundary profile instead of switching an initialized method.');
            end
            spec=method_sampler_profile(family,name);
        end
        function obj=apply_method_profile(obj,spec)
            if closed_enabled(obj) && ~strcmp(spec.name,'closed_ties_unit_fallback_v1')
                error('SuSClosedFallback:FreshModelRequired','Create a fresh model to change the closed-boundary method policy.');
            end
            if strcmp(spec.name,'closed_ties_unit_fallback_v1')
                if ~isfield(obj.G,'sus_runtime'),obj=closed_initialize(obj);
                elseif obj.G.sus_runtime.run_started
                    obj.G.sampler.aCS.Sigma0Mode=obj.G.sus_runtime.effective_sigma0_mode;
                end
            end
        end
        function obj = SuS_Bay_Nataf(LKF, method_profile)
            if nargin < 2 || isempty(method_profile), method_profile='author_latest'; end
            % Constructor
            % INPUT:
            %   LKF : likelihood function object
            obj.LKF = LKF;
            obj.Ndim = LKF.Ndim;
            obj=semcs.configure_method_profile(obj,'sus_bayes',method_profile);
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
            if closed_enabled(obj)
                Li = log(double(closed_contains(L,L_i)));
            elseif isinf(L_i)
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
            if closed_enabled(obj)
                closed_validate_values(L);
                p_i1 = p_i + log(sum(closed_contains(L,L_i1)) / Nsam);
            else
                p_i1 = p_i + log(sum(L > L_i1) / Nsam);
            end

            z_i1 = p_i + logmean(min(logminus(L,L_i), logminus(L_i1,L_i)), "all");
            z_tmp1 = p_i + logmean(logminus(L,L_i), "all");
            z_sum1 = logplus(logplus(logminus(z_sum,z_tmp), z_i), z_tmp1);

            g = [L_i, p_i, L_i1, p_i1, z_i1, z_tmp1, z_sum1];
            obj.G.parameters = g;
            g = obj.G;
        end

        function [u_sed, y_sed] = SltSed(obj, u, y, Nsam)
            % Select seed samples for next iteration
            if closed_enabled(obj)
                [u_sed,y_sed]=closed_project(obj,u,y,Nsam);
                return
            end
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
        function [obj,Ncal]=after_advance(obj,x,y,g)
            Ncal=0;
            if closed_enabled(obj)
                obj=closed_plan(obj,x{end},y{end},obj.G.run.NumSam,numel(x));
            end
        end
        function [y,g,obj,Ncal]=advance(obj,x,y,p)
            [g,obj]=obj.UpdObj(y,p);Ncal=0;
        end
        function [xs,ys,steps]=select_seeds(obj,x,y,N)
            [xs,ys]=obj.SltSed(x,y,N);steps=round(N/size(xs,2));
        end
        function yes=stop(obj),yes=obj.FlgCvg;end

    end

    methods(Hidden)
        function yes=closed_policy_enabled(obj),yes=closed_enabled(obj);end
        function obj=prepare_closed_run(obj),obj=closed_prepare_run(obj);end
    end

end

function yes=closed_enabled(b)
    yes=isfield(b.G,'boundary_rule')&&strcmp(b.G.boundary_rule,'closed')&& ...
        isfield(b.G,'seed_scale_policy')&&strcmp(b.G.seed_scale_policy,'unit_if_all_selected_reference_vectors_identical_else_baseline');
end
function b=closed_initialize(b)
    b.G.boundary_rule='closed';
    b.G.seed_scale_policy='unit_if_all_selected_reference_vectors_identical_else_baseline';
    b.G.baseline_sigma0_mode=[];
    b.G.sus_runtime=struct('run_started',false,'sampled_levels',0,'chain_plan',[], ...
        'last_seed_selection',[],'degenerate_selected_seed_set',false,'effective_sigma0_mode',[]);
end
function b=closed_prepare_run(b)
    assert(closed_enabled(b),'SuSClosedFallback:Profile','Explicit opt-in profile required.');
    assert(~b.G.sus_runtime.run_started&&b.G.sus_runtime.sampled_levels==0, ...
        'SuSClosedFallback:FreshModelRequired','Run this variant from a fresh initialized model.');
    validateattributes(b.G.run.NumSam,{'numeric'},{'scalar','positive','integer','finite'});
    assert(isempty(b.G.run.schedule),'SuSClosedFallback:Schedule','This minimum variant does not support prescribed G schedules.');
    assert(isequal(b.G.parameters,[-inf,0,-inf,0,-inf,-inf,-inf]), ...
        'SuSClosedFallback:FreshModelRequired','Initial parameters must be the ordinary fresh prior state.');
    base=b.G.sampler.aCS.Sigma0Mode;
    assert(any(strcmp(base,{'seed_std','unit'})),'SuSClosedFallback:Scale','Choose seed_std or unit.');
    b.G.baseline_sigma0_mode=base;
    b.G.sus_runtime.effective_sigma0_mode=base;b.G.sus_runtime.run_started=true;
end
function mask=closed_contains(logL,cut)
    if isinf(cut)&&cut<0,mask=true(size(logL));else,mask=logL>=cut;end
end
function closed_validate_values(logL)
    assert(isreal(logL)&&~isempty(logL)&&all(~isnan(logL(:)))&&all(logL(:)~=inf)&&any(isfinite(logL(:))), ...
        'SuSClosedFallback:Likelihood','Cached log likelihoods must exclude NaN/+Inf and contain a finite value.');
end
function z=closed_reference(b,u)
    opts=semcs.kernel_options(b,'aCS');
    if strcmp(opts.CoordinateMode,'model')
        z=b.to_reference(u);
    else
        L=eye(b.Ndim);
        if isfield(b.G,'L'),L=b.G.L;elseif isfield(b.G,'L_pi'),L=b.G.L_pi;end
        if isscalar(L),L=L*eye(b.Ndim);end
        assert(isequal(size(L),[b.Ndim,b.Ndim])&&all(isfinite(L(:)))&&rank(L)==b.Ndim, ...
            'SuSClosedFallback:Reference','Gaussian reference factor must be finite and nonsingular.');
        z=L\u;
    end
    assert(all(isfinite(z(:))),'SuSClosedFallback:Reference','Selected reference seed vectors must be finite.');
end
function b=closed_plan(b,u,y,N,stored_level)
    assert(b.G.sus_runtime.run_started,'SuSClosedFallback:Start','Capture the actual sampler baseline before selecting seeds.');
    validateattributes(N,{'numeric'},{'scalar','positive','integer','finite'});
    u=u(1:b.Ndim,:);y=y(1:b.Nfun,:);
    L=y(3,:);closed_validate_values(L);
    cut=b.G.parameters(3);eligible=find(closed_contains(L,cut));
    counts=round(N./round(N./(1:N)));
    available=find(counts<=numel(eligible),1,'last');
    if isempty(available),nc=0;chosen=[];positions=[];ns=inf;
    else
        nc=counts(available);positions=randperm(numel(eligible),nc);
        chosen=eligible(positions);ns=round(N/nc);
    end
    allsame=false;
    if nc>0
        z=closed_reference(b,u(:,chosen));
        allsame=all(z==z(:,1),'all');
    end
    base=b.G.baseline_sigma0_mode;effective=base;
    if strcmp(base,'seed_std')&&allsame,effective='unit';end
    b.G.sampler.aCS.Sigma0Mode=effective;
    event=struct('index_base',1,'eligible_indices',eligible,'selected_positions',positions, ...
        'selected_indices',chosen,'from_stored_level',stored_level,'source_count',numel(L), ...
        'target_count',N,'closed_hit_fraction',numel(eligible)/numel(L), ...
        'selected_reference_vectors_all_identical',allsame,'baseline_sigma0_mode',base, ...
        'effective_sigma0_mode',effective,'log_likelihood_threshold',cut);
    cp=struct('Nc',nc,'Ns',ns,'actual_count',nc*ns,'target_count',N, ...
        'source_stored_level',stored_level,'source_count',numel(L),'selected_indices',chosen, ...
        'parameters',b.G.parameters,'source_digest',closed_digest(u,y));
    if nc==0,cp.actual_count=0;end
    b.G.sus_runtime.chain_plan=cp;b.G.sus_runtime.last_seed_selection=event;
    b.G.sus_runtime.sampled_levels=stored_level;
    b.G.sus_runtime.degenerate_selected_seed_set=allsame;
    b.G.sus_runtime.effective_sigma0_mode=effective;
end
function [xs,ys]=closed_project(b,u,y,N)
    rt=b.G.sus_runtime;cp=rt.chain_plan;
    assert(~isempty(cp)&&cp.target_count==N&&isequaln(cp.parameters,b.G.parameters), ...
        'SuSClosedFallback:StalePlan','Prepare a seed plan under the current G and budget first.');
    u=u(1:b.Ndim,:);y=y(1:b.Nfun,:);
    assert(strcmp(cp.source_digest,closed_digest(u,y)), ...
        'SuSClosedFallback:StalePlan','The seed plan belongs to a different cached sample batch.');
    assert(strcmp(b.G.sampler.aCS.Sigma0Mode,rt.effective_sigma0_mode), ...
        'SuSClosedFallback:StalePlan','The sampler scale changed after its seed plan was prepared.');
    xs=u(:,cp.selected_indices);ys=y(:,cp.selected_indices);
    if ~isempty(ys),ys=b.UpdY(ys);end
    assert(isempty(ys)||all(isfinite(ys(1,:))), ...
        'SuSClosedFallback:SeedSupport','Every refreshed seed must belong to the current closed support.');
end
function value=closed_digest(varargin)
    d=java.security.MessageDigest.getInstance('SHA-256');
    for k=1:nargin
        a=varargin{k};
        assert(isa(a,'double')&&isreal(a),'SuSClosedFallback:Cache','Cached batches must be real double arrays.');
        shape=uint64([ndims(a),size(a)]);d.update(typecast(shape(:),'uint8'));
        d.update(typecast(a(:),'uint8'));
    end
    value=lower(reshape(dec2hex(typecast(d.digest(),'uint8'),2).',1,[]));
end

function spec = method_sampler_profile(family, name)
% Method-owned sampling settings; no budgets, thresholds or stops are changed.
% 'author_latest' follows the user's selected latest local author source policy.
% It does not authenticate a historical table-generation commit.
if nargin < 2 || isempty(name), name='author_latest'; end
family=char(family); name=char(name);
if strcmp(name,'closed_ties_unit_fallback_v1')
    spec=method_sampler_profile(family,'author_latest');
    spec.name=name;
    spec.evidence.closed_boundary_variant=struct('status','explicit_robustness_variant_not_literal_strict_paper', ...
        'detail','Closed support and actual closed mass; unit only for an exactly constant selected reference seed set with a seed_std baseline. No completion or efficiency guarantee.');
    return
end
if strcmp(name,'paper_settings'), name='paper'; end
if ~any(strcmp(name,{'author_latest','paper','public_20261102'}))
    error('method_sampler_profile:Name','Choose author_latest, paper, paper_settings or public_20261102.');
end
if ~any(strcmp(family,{'sus_bayes'}))
    error('method_sampler_profile:Family','Unsupported method family for this toolbox.');
end
options=struct();
options.aCS=struct('Sigma0Mode','unit','CarryLambda',true);
options.ES=struct('ESMode','coordinate','Ntry',20);
default_sampler='aCS';
support='public20261102 unit/carry aCS and coordinate ES replay';
if any(strcmp(family,{'gss_mis'})), default_sampler='CW_aMH'; end
if strcmp(family,'semis_bayes'), default_sampler='ES'; end
if strcmp(family,'resus'), default_sampler='ES'; end
if strcmp(name,'paper')
    switch family
        case 'semis_rel'
            options.aCS=struct('Sigma0Mode','unit','CarryLambda',false);
            support='retained SeMIS reliability src v1 aCS: unit sigma0, lambda reset';
        case 'sus_bayes'
            options.aCS=struct('Sigma0Mode','seed_std','CarryLambda',false);
            support='CS-MH article family supplemented by retained author aCS seed std/reset';
        case 'semis_bayes'
            options.aCS=struct('Sigma0Mode','seed_std','CarryLambda',false);
            options.ES=struct('ESMode','global','Ntry',[]);
            support='paper global vector-angle ESS, until accepted; optional aCS author-source supplement';
        case 'resus'
            default_sampler='aCS';
            options.aCS=struct('Sigma0Mode','seed_std','CarryLambda',false);
            options.ES=struct('ESMode','global','Ntry',[]);
            support='paper permits seed std or unit; retained author aCS seed std/reset; paper global ESS';
        case 'gss_mis'
            support='retained author CW_aMH seed std; no lambda setting applies to CW_aMH';
        case 'sus_rel'
            support='optional scalar SuS keeps public aCS settings; no method-paper kernel certification';
    end
end
if strcmp(name,'author_latest')
    switch family
        case 'semis_rel'
            options.aCS=struct('Sigma0Mode','unit','CarryLambda',false);
            support='latest reliability src v1: unit aCS, lambda reset';
        case 'sus_bayes'
            options.aCS=struct('Sigma0Mode','seed_std','CarryLambda',false);
            support='latest Bayesian SuS author aCS: seed std, lambda reset';
        case 'semis_bayes'
            options.aCS=struct('Sigma0Mode','seed_std','CarryLambda',false);
            options.ES=struct('ESMode','coordinate','Ntry',20);
            support='latest SeMIS Bayesian author ES: coordinate, Ntry20; optional aCS seed std/reset';
        case 'resus'
            default_sampler='ES';
            options.aCS=struct('Sigma0Mode','seed_std','CarryLambda',false);
            options.ES=struct('ESMode','coordinate','Ntry',20);
            support='latest ReSuS author driver: ES coordinate/Ntry20 and hybrid midpoint-strict threshold; optional aCS seed std/reset';
        case 'gss_mis'
            support='latest coherent src_MIS: CW seed std; active-history balance MIS';
        case 'sus_rel'
            support='optional scalar SuS retains public aCS settings; outside certified latest SeMIS method scope';
    end
end
options.aCS.AcceptArithmetic='log';options.aCS.CoordinateMode='identity';
options.aCS.GroupMode='available';options.aCS.AcceptanceRateMode='attempts';options.aCS.AdaptEmptyGroups=false;
options.ES.AcceptArithmetic='log';options.ES.CoordinateMode='identity';options.ES.SeedUpdateMode='each_step';
options.CW_aMH=struct('AcceptArithmetic','log','PriorMode','gaussian','CoordinateMode','identity');
if strcmp(name,'author_latest')
    options.aCS.AcceptArithmetic='exp';options.ES.AcceptArithmetic='exp';options.CW_aMH=struct('AcceptArithmetic','exp');
    options.aCS.GroupMode='fixed';options.aCS.AcceptanceRateMode='attempts';options.aCS.AdaptEmptyGroups=true;
    if strcmp(family,'semis_rel')
        options.aCS.GroupMode='available';options.aCS.AcceptanceRateMode='chain_steps';options.aCS.AdaptEmptyGroups=false;
    end
    if strcmp(family,'resus')
        options.aCS.CoordinateMode='model';options.ES.CoordinateMode='model';options.ES.SeedUpdateMode='all_accepted';
    end
end
evidence=struct('support',support, ...
    'historical_table_generator_commit_authenticated',false, ...
    'unknown_paper_budget_or_tolerance_not_filled',true, ...
    'reference_factor_separate_from_kernel_base_scale',true, ...
    'author_source_arithmetic_not_asserted_identical',true);
% Field-level evidence describes support, not historical table authentication.
recorded=struct('status','public_20261102_recorded_policy', ...
    'detail','Retained public defaults; not automatically an article setting');
author=struct('status','author_source_not_explicit_in_paper', ...
    'detail','Retained author aCS source: Lambda .6, NumIte 5, SedSav false, AcpOpt .44; exact table-generation commit unverified');
notcertified=struct('status','optional_kernel_not_method_paper_certified', ...
    'detail','Optional retained kernel policy, outside the method-paper benchmark kernel');
na=struct('status','not_applicable','detail','No lambda adaptation/target/groups in CW_aMH');
evidence.aCS=struct('Sigma0Mode',recorded,'CarryLambda',recorded, ...
    'Lambda',author,'NumIte',author,'SedSav',author,'AcpOpt',author);
evidence.ES=struct('ESMode',recorded,'Ntry',recorded);
evidence.CW_aMH=struct('Sigma0Mode',struct('status','author_source', ...
    'detail','Retained CW_aMH recomputes seed standard deviation per sampling call'), ...
    'Lambda',na,'CarryLambda',na,'AcpOpt',na,'NumIte',na);
if strcmp(name,'paper')
    switch family
        case 'semis_rel'
            evidence.aCS.Sigma0Mode=struct('status','author_source_not_explicit_in_paper', ...
                'source','SeMIS reliability src v1/Support-Fun/aCS.m L58', ...
                'detail','Unit proposal base scale; reference covariance inflation is separate');
            evidence.aCS.CarryLambda=struct('status','author_source_not_explicit_in_paper', ...
                'source','SeMIS reliability src v1/Support-Fun/aCS.m L30,L103', ...
                'detail','Lambda starts from object option at each call; no adapted value written back');
        case 'sus_bayes'
            evidence.aCS.Sigma0Mode=struct('status','article_family_plus_author_source_supplement', ...
                'detail','Seed standard deviation in whitened reference coordinates; exact table-generating implementation unverified');
            evidence.aCS.CarryLambda=struct('status','author_source_not_explicit_in_paper', ...
                'detail','Reset from configured Lambda each call; no cross-level carry claimed from paper text');
        case 'semis_bayes'
            evidence.aCS.Sigma0Mode=notcertified;evidence.aCS.CarryLambda=notcertified;
        case 'resus'
            evidence.aCS.Sigma0Mode=struct('status','paper_allows_choice_author_source_selects', ...
                'source','Li et al. 2025 p3 Sec2.2.1; retained src/github/MCMC-samplers/aCS.m L76', ...
                'detail','Paper permits seed std or unit; this profile selects author-source seed std');
            evidence.aCS.CarryLambda=struct('status','author_source_not_explicit_in_paper', ...
                'source','retained src/github/MCMC-samplers/aCS.m L55,L113', ...
                'detail','Local lambda resets each call; final adaptation not returned to object');
            evidence.aCS.AcpOpt=struct('status','paper_explicit', ...
                'source','Li et al. 2025 p3 Sec2.2.1','detail','Target acceptance approximately .44');
        otherwise
            evidence.aCS.Sigma0Mode=notcertified;evidence.aCS.CarryLambda=notcertified;
    end
    if any(strcmp(family,{'semis_bayes','resus'}))
        if strcmp(family,'resus'),source='Li et al. 2025 p4 Eq12';
        else,source='SeMIS Bayesian paper global elliptical slice sampling description';end
        evidence.ES.ESMode=struct('status','paper_explicit','source',source, ...
            'detail','One angle for full vector, full 2pi bracket, shrink until accepted');
        evidence.ES.Ntry=struct('status','not_applicable','detail','Global paper ESS has no finite rejection-attempt cap');
    end
end
if strcmp(name,'author_latest')
    evidence.policy=struct('status','user_selected_latest_author_source', ...
        'detail','Selected by runnable entry/algorithm file mtime; not historical table identity');
    if any(strcmp(family,{'semis_rel','sus_bayes','semis_bayes','resus'}))
        evidence.aCS.Sigma0Mode=struct('status','latest_author_source', ...
            'detail',['Configured base scale: ' options.aCS.Sigma0Mode]);
        evidence.aCS.CarryLambda=struct('status','latest_author_source', ...
            'detail','Per-call lambda reset from configured option; final runtime value remains recorded');
    end
    if any(strcmp(family,{'semis_bayes','resus'}))
        evidence.ES.ESMode=struct('status','latest_author_source', ...
            'detail','Author coordinate-angle ES is retained instead of replacing it with a paper-formula policy');
        evidence.ES.Ntry=struct('status','latest_author_source','detail','Finite Ntry=20');
    end
end
if any(strcmp(name,{'paper','author_latest'}))&&strcmp(family,'sus_bayes')
    explicit=struct('status','paper_explicit', ...
        'source','SuS-Bay proof p6 Algorithm 2', ...
        'detail','Seed standard deviation; initial lambda .6; target acceptance .44. Exact five groups is an author-source choice within the printed adaptation-fraction range.');
    evidence.aCS.Sigma0Mode=explicit;evidence.aCS.Lambda=explicit;evidence.aCS.AcpOpt=explicit;
end
evidence.author_group_and_acceptance_arithmetic=struct('status','not_replayed_in_this_increment', ...
    'detail','author_latest selects retained source group/exp arithmetic; paper/public retain their recorded numerical policies');
spec=struct('family',family,'name',name,'default_sampler',default_sampler, ...
    'options',options,'evidence',evidence);
end
