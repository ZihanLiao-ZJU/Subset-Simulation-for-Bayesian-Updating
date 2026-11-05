% Adapted from Zihan Liao's public MATLAB MCMC samplers; repository LICENSE applies.
classdef CW_aMH < MCMCBase
    % Componentwise uniform random walk MH in whitened Gaussian coordinates.
    properties(Dependent)
        PriorMode
    end
    methods
        function value=get.PriorMode(obj),value=obj.read_option('PriorMode');end
        function obj=set.PriorMode(obj,value),obj=obj.write_option('PriorMode',value);end
        function obj=CW_aMH(b)
            if nargin<1, b=[]; end
            obj@MCMCBase(b);
        end
        function [u,y,Ncal,b]=SamGen(obj,varargin)
            [b,u0,y0,Ncal,save,L,Nc,Ns]=obj.prepare(varargin{:});opts=MCMCBase.kernel_options(b,'CW_aMH');
            d=b.Ndim; z0=L\u0; u=zeros(d,Nc,Ns); y=zeros(b.Nfun,Nc,Ns);
            du=2*rand(d,Nc,Ns)-1; p=rand(d,Nc,Ns); sigma=max(std(z0,0,2),eps);
            before=MCMCBase.profile_snapshot(b,'CW_aMH',struct('Sigma0Mode','seed_std','sigma0',sigma, ...
                'coordinates','whitened','LambdaApplicable',false,'AcceptanceTargetApplicable',false,'GroupsApplicable',false));
            for k=1:Ns
                if k==1 && save, u(:,:,k)=u0; y(:,:,k)=y0; continue; end
                z=z0+sigma.*du(:,:,k);
                if strcmp(opts.PriorMode,'model')
                    fresh=zeros(size(z));old=fresh;
                    for j=1:d,fresh(j,:)=b.marginal_logpdf(z(j,:),j);old(j,:)=b.marginal_logpdf(z0(j,:),j);end
                    rej=obj.reject(fresh,old,p(:,:,k));
                else,rej=obj.reject(-.5*z.^2,-.5*z0.^2,p(:,:,k));end
                z(rej)=z0(rej);
                ids=find(any(z~=z0,1));
                if ~isempty(ids)
                    proposal=L*z(:,ids); [b,yp,nn]=obj.evaluate(b,proposal); Ncal=Ncal+nn;
                    rej=obj.reject(yp(1,:),y0(1,ids),rand(1,numel(ids))); yes=~rej;
                    z0(:,ids(yes))=z(:,ids(yes)); u0(:,ids(yes))=proposal(:,yes); y0(:,ids(yes))=yp(:,yes);
                end
                u(:,:,k)=u0; y(:,:,k)=y0;
            end
            b.G.sampler_record=struct('kernel','CW_aMH','before',before,'after',before);
        end
    end
end
