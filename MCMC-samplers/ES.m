% Adapted from Zihan Liao's public MATLAB MCMC samplers; repository LICENSE applies.
classdef ES < MCMCBase
    % User's componentwise [0,pi] elliptical slice kernel, whitened for L.
    properties
        Ntry = 20
    end
    methods
        function obj=ES(b)
            if nargin<1, b=[]; end
            obj@MCMCBase(b);
        end
        function [u,y,Ncal,b]=SamGen(obj,varargin)
            [b,u0,y0,Ncal,save,L,Nc,Ns]=obj.prepare(varargin{:});
            d=b.Ndim; z0=L\u0; u=zeros(d,Nc,Ns); y=zeros(b.Nfun,Nc,Ns);
            v=randn(d,Nc,Ns); p=rand(1,Nc,Ns);
            for k=1:Ns
                if k==1 && save, u(:,:,k)=u0; y(:,:,k)=y0; continue; end
                R=sqrt(z0.^2+v(:,:,k).^2); ratio=zeros(size(z0));
                nz=R>0; ratio(nz)=z0(nz)./R(nz); oldangle=acos(max(-1,min(1,ratio)));
                ids=1:Nc; lo=zeros(d,Nc); hi=pi*ones(d,Nc); z=z0;
                u(:,:,k)=u0; y(:,:,k)=y0;
                for j=1:obj.Ntry
                    angle=lo+(hi-lo).*rand(d,numel(ids)); proposed=cos(angle).*R(:,ids);
                    proposal=L*proposed; [b,yp,nn]=obj.evaluate(b,proposal); Ncal=Ncal+nn;
                    rej=obj.reject(yp(1,:),y0(1,ids),p(1,ids,k)); yes=~rej;
                    z(:,ids(yes))=proposed(:,yes); u(:,ids(yes),k)=proposal(:,yes); y(:,ids(yes),k)=yp(:,yes);
                    if ~any(rej), break; end
                    ids=ids(rej); angle=angle(:,rej); lo=lo(:,rej); hi=hi(:,rej);
                    upper=angle>oldangle(:,ids); lo(~upper)=angle(~upper); hi(upper)=angle(upper);
                end
                z0=z; u0=u(:,:,k); y0=y(:,:,k);
            end
        end
    end
end
