% Adapted from Zihan Liao's public MATLAB MCMC samplers; repository LICENSE applies.
classdef CW_aMH < MCMCBase
    % Componentwise uniform random walk MH in whitened Gaussian coordinates.
    methods
        function obj=CW_aMH(b)
            if nargin<1, b=[]; end
            obj@MCMCBase(b);
        end
        function [u,y,Ncal,b]=SamGen(obj,varargin)
            [b,u0,y0,Ncal,save,L,Nc,Ns]=obj.prepare(varargin{:});
            d=b.Ndim; z0=L\u0; u=zeros(d,Nc,Ns); y=zeros(b.Nfun,Nc,Ns);
            du=2*rand(d,Nc,Ns)-1; p=rand(d,Nc,Ns); sigma=max(std(z0,0,2),eps);
            for k=1:Ns
                if k==1 && save, u(:,:,k)=u0; y(:,:,k)=y0; continue; end
                z=z0+sigma.*du(:,:,k); rej=obj.reject(-.5*z.^2,-.5*z0.^2,p(:,:,k)); z(rej)=z0(rej);
                ids=find(any(z~=z0,1));
                if ~isempty(ids)
                    proposal=L*z(:,ids); [b,yp,nn]=obj.evaluate(b,proposal); Ncal=Ncal+nn;
                    rej=obj.reject(yp(1,:),y0(1,ids),rand(1,numel(ids))); yes=~rej;
                    z0(:,ids(yes))=z(:,ids(yes)); u0(:,ids(yes))=proposal(:,yes); y0(:,ids(yes))=yp(:,yes);
                end
                u(:,:,k)=u0; y(:,:,k)=y0;
            end
        end
    end
end
