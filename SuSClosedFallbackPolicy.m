classdef SuSClosedFallbackPolicy
    % Method-owned closed support and next-seed plan. No sampling kernel here.
    methods(Static)
        function yes=enabled(b)
            yes=isfield(b.G,'boundary_rule')&&strcmp(b.G.boundary_rule,'closed')&& ...
                isfield(b.G,'seed_scale_policy')&&strcmp(b.G.seed_scale_policy,'unit_if_all_selected_reference_vectors_identical_else_baseline');
        end
        function b=initialize(b)
            b.G.boundary_rule='closed';
            b.G.seed_scale_policy='unit_if_all_selected_reference_vectors_identical_else_baseline';
            b.G.baseline_sigma0_mode=[];
            b.G.sus_runtime=struct('run_started',false,'sampled_levels',0,'chain_plan',[], ...
                'last_seed_selection',[],'degenerate_selected_seed_set',false,'effective_sigma0_mode',[]);
        end
        function b=prepare_run(b)
            assert(SuSClosedFallbackPolicy.enabled(b),'SuSClosedFallback:Profile','Explicit opt-in profile required.');
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
        function mask=contains(logL,cut)
            if isinf(cut)&&cut<0,mask=true(size(logL));else,mask=logL>=cut;end
        end
        function validate_values(logL)
            assert(isreal(logL)&&~isempty(logL)&&all(~isnan(logL(:)))&&all(logL(:)~=inf)&&any(isfinite(logL(:))), ...
                'SuSClosedFallback:Likelihood','Cached log likelihoods must exclude NaN/+Inf and contain a finite value.');
        end
        function z=reference(b,u)
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
        function b=plan(b,u,y,N,stored_level)
            assert(b.G.sus_runtime.run_started,'SuSClosedFallback:Start','Capture the actual sampler baseline before selecting seeds.');
            validateattributes(N,{'numeric'},{'scalar','positive','integer','finite'});
            u=u(1:b.Ndim,:);y=y(1:b.Nfun,:);
            L=y(3,:);SuSClosedFallbackPolicy.validate_values(L);
            cut=b.G.parameters(3);eligible=find(SuSClosedFallbackPolicy.contains(L,cut));
            counts=round(N./round(N./(1:N)));
            available=find(counts<=numel(eligible),1,'last');
            if isempty(available),nc=0;chosen=[];positions=[];ns=inf;
            else
                nc=counts(available);positions=randperm(numel(eligible),nc);
                chosen=eligible(positions);ns=round(N/nc);
            end
            allsame=false;
            if nc>0
                z=SuSClosedFallbackPolicy.reference(b,u(:,chosen));
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
                'parameters',b.G.parameters,'source_digest',SuSClosedFallbackPolicy.digest(u,y));
            if nc==0,cp.actual_count=0;end
            b.G.sus_runtime.chain_plan=cp;b.G.sus_runtime.last_seed_selection=event;
            b.G.sus_runtime.sampled_levels=stored_level;
            b.G.sus_runtime.degenerate_selected_seed_set=allsame;
            b.G.sus_runtime.effective_sigma0_mode=effective;
        end
        function [xs,ys]=project(b,u,y,N)
            rt=b.G.sus_runtime;cp=rt.chain_plan;
            assert(~isempty(cp)&&cp.target_count==N&&isequaln(cp.parameters,b.G.parameters), ...
                'SuSClosedFallback:StalePlan','Prepare a seed plan under the current G and budget first.');
            u=u(1:b.Ndim,:);y=y(1:b.Nfun,:);
            assert(strcmp(cp.source_digest,SuSClosedFallbackPolicy.digest(u,y)), ...
                'SuSClosedFallback:StalePlan','The seed plan belongs to a different cached sample batch.');
            assert(strcmp(b.G.sampler.aCS.Sigma0Mode,rt.effective_sigma0_mode), ...
                'SuSClosedFallback:StalePlan','The sampler scale changed after its seed plan was prepared.');
            xs=u(:,cp.selected_indices);ys=y(:,cp.selected_indices);
            if ~isempty(ys),ys=b.UpdY(ys);end
            assert(isempty(ys)||all(isfinite(ys(1,:))), ...
                'SuSClosedFallback:SeedSupport','Every refreshed seed must belong to the current closed support.');
        end
        function value=digest(varargin)
            d=java.security.MessageDigest.getInstance('SHA-256');
            for k=1:nargin
                a=varargin{k};
                assert(isa(a,'double')&&isreal(a),'SuSClosedFallback:Cache','Cached batches must be real double arrays.');
                shape=uint64([ndims(a),size(a)]);d.update(typecast(shape(:),'uint8'));
                d.update(typecast(a(:),'uint8'));
            end
            value=lower(reshape(dec2hex(typecast(d.digest(),'uint8'),2).',1,[]));
        end
    end
end
