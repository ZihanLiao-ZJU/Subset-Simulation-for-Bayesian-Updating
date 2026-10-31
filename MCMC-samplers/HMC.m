% Adapted from Zihan Liao's public MATLAB MCMC samplers; repository LICENSE applies.
classdef HMC < MCMCBase
    % Exact Gaussian-reference HMC rotation with likelihood MH correction.
    % Rejection-based Gaussian HMC, not arbitrary-potential gradient HMC.
    properties
        InitialTime = pi/4
        NumIte = 5
        Adapt = true
    end
    methods
        function obj=HMC(b)
            if nargin<1, b=[]; end
            obj@MCMCBase(b);
        end
        function [u,y,Ncal,b]=SamGen(obj,varargin)
            [b,u0,y0,Ncal,save,L,Nc,Ns]=obj.prepare(varargin{:});
            d=b.Ndim; u=zeros(d,Nc,Ns); y=zeros(b.Nfun,Nc,Ns);
            v=reshape(L*reshape(randn(d,Nc,Ns),d,[]),d,Nc,Ns); p=rand(1,Nc,Ns);
            time=obj.state(b,'HMC_time',obj.InitialTime);
            if ~isscalar(time)||~isfinite(time)||time<=0||time>pi/2, error('HMC:State','G.HMC_time must lie in (0,pi/2].'); end
            width=max(ceil(Nc/obj.NumIte),1);
            for start=1:width:Nc
                ids=start:min(start+width-1,Nc); count=0; attempts=0;
                for k=1:Ns
                    if k==1 && save, u(:,ids,k)=u0(:,ids); y(:,ids,k)=y0(:,ids); continue; end
                    proposal=cos(time)*u0(:,ids)+sin(time)*v(:,ids,k);
                    [b,yp,nn]=obj.evaluate(b,proposal); Ncal=Ncal+nn;
                    rej=obj.reject(yp(1,:),y0(1,ids),p(1,ids,k)); yes=~rej;
                    u0(:,ids(yes))=proposal(:,yes); y0(:,ids(yes))=yp(:,yes);
                    u(:,ids,k)=u0(:,ids); y(:,ids,k)=y0(:,ids);
                    count=count+sum(yes); attempts=attempts+numel(ids);
                end
                if attempts>0 && obj.Adapt
                    rate=count/attempts;
                    if rate<.3, time=asin(min(1,sin(time)*exp((rate-.3)/2)));
                    elseif rate>.5, time=asin(min(1,sin(time)*exp((rate-.5)/2))); end
                end
            end
            b.G.HMC_time=time;
        end
    end
end
