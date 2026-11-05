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
            b=semcs.configure(b);
            kernel=semcs.kernel_name(obj);
            for key=fieldnames(obj.UnboundOptions)'
                b.G.sampler.(kernel).(key{1})=obj.UnboundOptions.(key{1});
            end
            obj.ModelState=b;obj.UnboundOptions=struct();
        end
        function opts=configuration_options(obj)
            if isempty(obj.ModelState),opts=obj.UnboundOptions;
            else,opts=obj.ModelState.G.sampler.(semcs.kernel_name(obj));end
        end
        function obj=unbound_prototype(obj)
            % Engine stores kernel identity only, never a second model G.
            obj.ModelState=[];obj.UnboundOptions=struct();
        end
        function value=read_option(obj,name)
            defaults=semcs.sampler_defaults(semcs.kernel_name(obj));value=defaults.(name);
            opts=obj.configuration_options();
            if isfield(opts,name),value=opts.(name);end
        end
        function obj=write_option(obj,name,value)
            if isempty(obj.ModelState),obj.UnboundOptions.(name)=value;
            else,obj.ModelState.G=semcs.setfield(obj.ModelState.G,['sampler.' semcs.kernel_name(obj) '.' name],value);end
        end
        function d=get.Ndim(obj), d=obj.BayStc.Ndim; end
        function n=get.Nfun(obj), n=obj.BayStc.Nfun; end
        function obj=MCMCBase(b)
            if nargin>0, obj.BayStc=b; end
        end
        function [b,u,y,n,save,L,Nc,Ns]=prepare(obj,varargin)
            b=semcs.configure(obj.BayStc); d=b.Ndim; L=eye(d);
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
                u=varargin{1}; y=varargin{2}; Ns=varargin{3}; Nc=size(u,2); n=0; opts=semcs.kernel_options(b,semcs.kernel_name(obj));save=opts.SedSav;
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
            opts=semcs.kernel_options(obj.BayStc,semcs.kernel_name(obj));
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
        function value=option(obj,b,name),opts=semcs.kernel_options(b,semcs.kernel_name(obj));value=opts.(name);end

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

end
