% Adapted from Zihan Liao's public MATLAB MCMC samplers; repository LICENSE applies.
classdef aCS < MCMCBase
    % sigma0=1; updated lambda is returned in b.G.lambda across levels.
    properties
        Lambda = 0.6
        AcpOpt = 0.44
        NumIte = 5
        Adapt = true
    end
    methods
        function obj=aCS(b)
            if nargin<1, b=[]; end
            obj@MCMCBase(b);
        end
        function [u,y,Ncal,b]=SamGen(obj,varargin)
            [b,u0,y0,Ncal,save,L,Nc,Ns]=obj.prepare(varargin{:});
            d=b.Ndim; u=zeros(d,Nc,Ns); y=zeros(b.Nfun,Nc,Ns);
            v=reshape(L*reshape(randn(d,Nc,Ns),d,[]),d,Nc,Ns); p=rand(1,Nc,Ns);
            lambda=obj.state(b,'lambda',obj.Lambda);
            if ~isscalar(lambda)||~isfinite(lambda)||lambda<=0, error('aCS:State','G.lambda must be finite and positive.'); end
            width=max(ceil(Nc/obj.NumIte),1); group=0;
            for start=1:width:Nc
                group=group+1; ids=start:min(start+width-1,Nc);
                sigma=min(lambda,1); rho=sqrt(1-sigma^2); count=0; attempts=0;
                for k=1:Ns
                    if k==1 && save, u(:,ids,k)=u0(:,ids); y(:,ids,k)=y0(:,ids); continue; end
                    proposal=rho*u0(:,ids)+sigma*v(:,ids,k);
                    [b,yp,nn]=obj.evaluate(b,proposal); Ncal=Ncal+nn;
                    rej=obj.reject(yp(1,:),y0(1,ids),p(1,ids,k)); yes=~rej;
                    u0(:,ids(yes))=proposal(:,yes); y0(:,ids(yes))=yp(:,yes);
                    u(:,ids,k)=u0(:,ids); y(:,ids,k)=y0(:,ids);
                    count=count+sum(yes); attempts=attempts+numel(ids);
                end
                if attempts>0 && obj.Adapt, lambda=exp(log(lambda)+(count/attempts-obj.AcpOpt)/sqrt(group)); end
            end
            b.G.lambda=lambda;
        end
    end
end
