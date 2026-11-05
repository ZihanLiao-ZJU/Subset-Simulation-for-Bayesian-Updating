classdef MCMCBase
    properties(Dependent)
        SedSav
        Lambda
        AcpOpt
        NumIte
        Adapt
        Ntry
        InitialTime
        AcceptArithmetic
        CoordinateMode
    end

    % Common Gaussian/Nataf target protocol: prior times likelihood factor.
    properties(Dependent)
        BayStc
    end
    properties(Access=private)
        ModelState=[]
        UnboundOptions=struct()
    end
    properties (Dependent)
        Ndim
        Nfun
    end
    methods
        function b=get.BayStc(obj),b=obj.ModelState;end
        function obj=set.BayStc(obj,b)
            if isempty(b)
                if ~isempty(obj.ModelState),obj.UnboundOptions=obj.configuration_options();end
                obj.ModelState=[];
                return
            end
            b=MCMCBase.configure(b);
            kernel=MCMCBase.kernel_name(obj);
            for key=fieldnames(obj.UnboundOptions)'
                b.G.sampler.(kernel).(key{1})=obj.UnboundOptions.(key{1});
            end
            obj.ModelState=b;obj.UnboundOptions=struct();
        end
        function opts=configuration_options(obj)
            if isempty(obj.ModelState),opts=obj.UnboundOptions;
            else,opts=obj.ModelState.G.sampler.(MCMCBase.kernel_name(obj));end
        end
        function obj=unbound_prototype(obj)
            % Engine stores kernel identity only, never a second model G.
            obj.ModelState=[];obj.UnboundOptions=struct();
        end
        function value=read_option(obj,name)
            defaults=MCMCBase.sampler_defaults(MCMCBase.kernel_name(obj));value=defaults.(name);
            opts=obj.configuration_options();
            if isfield(opts,name),value=opts.(name);end
        end
        function obj=write_option(obj,name,value)
            if isempty(obj.ModelState),obj.UnboundOptions.(name)=value;
            else,obj.ModelState.G=MCMCBase.state_setfield(obj.ModelState.G,['sampler.' MCMCBase.kernel_name(obj) '.' name],value);end
        end
        function d=get.Ndim(obj), d=obj.BayStc.Ndim; end
        function n=get.Nfun(obj), n=obj.BayStc.Nfun; end
        function obj=MCMCBase(b)
            if nargin>0, obj.BayStc=b; end
        end
        function [b,u,y,n,save,L,Nc,Ns]=prepare(obj,varargin)
            b=MCMCBase.configure(obj.BayStc); d=b.Ndim; L=eye(d);
            if isfield(b.G,'L'), L=b.G.L;
            elseif isfield(b.G,'L_pi'), L=b.G.L_pi; end
            if isscalar(L), L=L*eye(d); end
            if ~isequal(size(L),[d,d]) || any(~isfinite(L(:))) || rank(L)~=d
                error('MCMC:Reference','Gaussian factor must be finite nonsingular Ndim by Ndim.');
            end
            if numel(varargin)==2
                Nc=varargin{1}; Ns=varargin{2};
                if ~isscalar(Nc) || ~isscalar(Ns) || ~isfinite(Nc) || ~isfinite(Ns) || Nc<1 || Ns<1 || Nc~=floor(Nc) || Ns~=floor(Ns), error('MCMC:Count','Counts must be positive integers.'); end
                if ismethod(b,'SamGen'), u=b.SamGen(Nc); else, u=L*randn(d,Nc); end
                [b,y,n]=obj.evaluate(b,u); save=true;
            elseif numel(varargin)==3
                u=varargin{1}; y=varargin{2}; Ns=varargin{3}; Nc=size(u,2); n=0; opts=MCMCBase.kernel_options(b,MCMCBase.kernel_name(obj));save=opts.SedSav;
            else
                error('MCMC:Arguments','SamGen(Nchains,Nsteps) or SamGen(U_seed,Y_seed,Nsteps).');
            end
            if Nc<1 || Ns<1 || Nc~=floor(Nc) || Ns~=floor(Ns) || ~isequal(size(u),[d,Nc]) || ~isequal(size(y),[b.Nfun,Nc])
                error('MCMC:Shape','Invalid chain count, step count or packed seed outputs.');
            end
        end
        function [b,y,n]=evaluate(~,b,u)
            b=b.EvlY(u); y=[b.EvlLKF;b.EvlPDF;b.Y]; n=b.Ncal;
            if ~isequal(size(y),[b.Nfun,size(u,2)])
                error('MCMC:Output','Invalid packed target outputs.');
            end
        end
        function r=reject(obj,new,old,p)
            opts=MCMCBase.kernel_options(obj.BayStc,MCMCBase.kernel_name(obj));
            if strcmp(opts.AcceptArithmetic,'exp'),r=exp(new-old)<p;
            else,delta=new-old;r=isnan(delta)|delta<log(max(p,realmin));end
        end
        function z=to_reference(~,b,u,L,opts)
            if strcmp(opts.CoordinateMode,'model'),z=b.to_reference(u);else,z=L\u;end
        end
        function u=from_reference(~,b,z,L,opts)
            if strcmp(opts.CoordinateMode,'model'),u=b.from_reference(z);else,u=L*z;end
        end
        function value=state(~,b,key,default)
            value=default; if isfield(b.G,key) && ~isempty(b.G.(key)), value=b.G.(key); end
        end
    end
    methods
        function value=get.AcceptArithmetic(obj),value=obj.read_option('AcceptArithmetic');end
        function obj=set.AcceptArithmetic(obj,value),obj=obj.write_option('AcceptArithmetic',value);end
        function value=get.CoordinateMode(obj),value=obj.read_option('CoordinateMode');end
        function obj=set.CoordinateMode(obj,value),obj=obj.write_option('CoordinateMode',value);end

    end

    methods
        function value=get.SedSav(obj)
            value=obj.read_option('SedSav');
        end
        function obj=set.SedSav(obj,value)
            obj=obj.write_option('SedSav',value);
        end
        function value=get.Lambda(obj)
            value=obj.read_option('Lambda');
        end
        function obj=set.Lambda(obj,value)
            obj=obj.write_option('Lambda',value);
        end
        function value=get.AcpOpt(obj)
            value=obj.read_option('AcpOpt');
        end
        function obj=set.AcpOpt(obj,value)
            obj=obj.write_option('AcpOpt',value);
        end
        function value=get.NumIte(obj)
            value=obj.read_option('NumIte');
        end
        function obj=set.NumIte(obj,value)
            obj=obj.write_option('NumIte',value);
        end
        function value=get.Adapt(obj)
            value=obj.read_option('Adapt');
        end
        function obj=set.Adapt(obj,value)
            obj=obj.write_option('Adapt',value);
        end
        function value=get.Ntry(obj)
            value=obj.read_option('Ntry');
        end
        function obj=set.Ntry(obj,value)
            obj=obj.write_option('Ntry',value);
        end
        function value=get.InitialTime(obj)
            value=obj.read_option('InitialTime');
        end
        function obj=set.InitialTime(obj,value)
            obj=obj.write_option('InitialTime',value);
        end

    end


    methods(Static,Hidden)
        function [b,prototype]=adopt_sampler(b,s)
        % Import explicitly supplied kernel options, keeping model/run/live G state.
        assert(~isa(s,'handle'),'semcs:HandleSamplerUnsupported', ...
            'Use a value-class sampler; a handle sampler cannot be detached without mutating its caller.');
        b=MCMCBase.configure(b);
        if isa(s,'MCMCBase')
            opts=s.configuration_options();
            kernel=MCMCBase.kernel_name(s);
            for key=fieldnames(opts)'
                b.G.sampler.(kernel).(key{1})=opts.(key{1});
            end
            prototype=s.unbound_prototype();
        else
            % Value-class custom samplers retain the writable BayStc/SamGen protocol.
            prototype=s;prototype.BayStc=[];
        end
        end
        function b=configure(b)
        defaults=struct('NumSam',1000,'MaxIte',10000,'Verbose',true, ...
         'stop_before_seeds',false,'select_at_cap',false,'level_offset',0,'schedule',{{}});
        if isprop(b,'Nsam'),defaults.NumSam=b.Nsam;end
        if ~isfield(b.G,'p'),b.G.p=0.1;end
        for key=fieldnames(defaults)'
         name=key{1};if ~isfield(b.G,'run')||~isfield(b.G.run,name),b.G=MCMCBase.state_setfield(b.G,['run.' name],defaults.(name));end
        end
        if ~isfield(b.G,'sampler'),b.G.sampler=struct();end
        for kernel={'aCS','ES','HMC','CW_aMH'}
         name=kernel{1};opts=MCMCBase.sampler_defaults(name);
         if isfield(b.G.sampler,name)
          incoming=b.G.sampler.(name);
          for key=fieldnames(incoming)',opts.(key{1})=incoming.(key{1});end
         end
         b.G.sampler.(name)=opts;
        end
        end
        function b=configure_method_profile(b,family,name)
        % Method-owned sampling configuration only: no threshold/budget/stop changes.
        if nargin<3,name=[];end
        assert(ismethod(b,'sampler_profile'),'semcs:Profile','The model must define sampler_profile(family,name).');
        b=MCMCBase.configure(b);spec=b.sampler_profile(family,name);
        assert(isstruct(spec)&&isscalar(spec)&&all(isfield(spec,{'family','name','default_sampler','options'})), ...
         'semcs:Profile','The method provider must return one complete profile spec.');
        for field={'family','name','default_sampler'}
         value=spec.(field{1});
         assert((ischar(value)&&isrow(value)&&~isempty(value))||(isstring(value)&&isscalar(value)&&strlength(value)>0), ...
          'semcs:Profile','Profile identifiers must be nonempty text.');
        end
        assert(isstruct(spec.options)&&isscalar(spec.options), ...
         'semcs:Profile','Profile options must be one scalar struct.');
        for kernel=fieldnames(spec.options)'
         key=kernel{1};patch=spec.options.(key);
         assert(isstruct(patch)&&isscalar(patch),'semcs:Profile','Each kernel options block must be one struct.');
         for option=fieldnames(patch)',b.G.sampler.(key).(option{1})=patch.(option{1});end
        end
        b.G.method_family=spec.family;b.G.method_profile=spec.name;
        b.G.default_sampler=spec.default_sampler;
        end
        function name=kernel_name(s)
        if isa(s,'aCS'),name='aCS';
        elseif isa(s,'ES'),name='ES';
        elseif isa(s,'HMC'),name='HMC';
        elseif isa(s,'CW_aMH'),name='CW_aMH';
        else,name=matlab.lang.makeValidName(class(s));end
        end
        function opts=kernel_options(b,kernel)
        opts=b.G.sampler.(kernel);
        end
        function opts=sampler_defaults(kernel)
        opts=struct('SedSav',false,'Lambda',0.6,'AcpOpt',0.44,'NumIte',5,'Adapt',true,'Ntry',20,'InitialTime',pi/4);
        opts.AcceptArithmetic='log';opts.CoordinateMode='identity';
        switch kernel
         case 'aCS',opts.Sigma0Mode='unit';opts.CarryLambda=true;opts.GroupMode='available';opts.AcceptanceRateMode='attempts';
         case 'ES',opts.ESMode='coordinate';opts.SeedUpdateMode='each_step';
         case 'CW_aMH',opts.Sigma0Mode='seed_std';opts.PriorMode='gaussian';
        end
        end
        function g=state_setfield(g,path,value)
        parts=strsplit(path,'.');g=MCMCBase.state_put(g,parts,value);
        end
        function g=state_put(g,p,value)
        if ~isstruct(g),g=struct();end
        if numel(p)==1,g.(p{1})=value;return;end
        if ~isfield(g,p{1}),g.(p{1})=struct();end
        g.(p{1})=MCMCBase.state_put(g.(p{1}),p(2:end),value);
        end
    end
end
