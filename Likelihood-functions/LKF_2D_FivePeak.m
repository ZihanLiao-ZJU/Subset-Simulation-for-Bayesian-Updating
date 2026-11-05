classdef LKF_2D_FivePeak
    % Paper Bayesian example B: five unnormalized Gaussian-shaped peaks.
    properties
        Ndim=2
        Pdis
    end
    methods
        function obj=LKF_2D_FivePeak()
            obj.Pdis=Nataf(repmat(makedist('Uniform','lower',-1,'upper',1),2,1),eye(2));
        end
        function value=EvlLKF(~,theta)
            mu=[-.4 -.4;-.35 .2;-.2 .15;.1 -.15;.45 .1];
            sigma=[.01 .01 .03 .02 .05];amplitude=[.5 1 .8 .5 .6];
            terms=zeros(5,size(theta,2));
            for j=1:5,terms(j,:)=log(amplitude(j))-sum((theta-mu(j,:)').^2,1)/(2*sigma(j)^2);end
            top=max(terms,[],1);value=top+log(sum(exp(terms-top),1));
        end
    end
end
