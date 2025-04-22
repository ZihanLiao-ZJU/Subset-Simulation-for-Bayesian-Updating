function [g,q] = ThdGen(IntL,IntP,y,p)
% Generate the threshods in subset simulation
% Syntax:
% ----------------------------------------------------------------------
% [g,q] = ThdGen(IntL,IntP,y,p)
% ----------------------------------------------------------------------
% INPUTS:
% IntL : Intermediate likelihood function
% IntP : Intermediate prior distribution
% y    : target function value of u                            [6,Nsam]
%       --intermediate likelihood function Li
%       --likelihood function L
%       --target function f
%       --intermediate prior PDF pi
%       --prior PDF p
%       --reference distribution PDF q
% p    : filtering ratio
% ----------------------------------------------------------------------
% OUTPUTS:
% g    : threshold for IntL
% q    : threshold for IntP
% ----------------------------------------------------------------------
Nsze = size(y);
Nsam = prod(Nsze(2:end));
G = IntL.G;
Q = IntP.Q;
err_tol = 1e-3;
maxi = 100;
options = optimset('Display', 'off');


if IntP.FlgThd && ~IntL.FlgThd
    % if there exsits a threshold generation function in IntP
    q = Intp.GenThd(y,p);
    g = G;
elseif ~IntP.FlgThd && IntL.FlgThd
    % if there exsits a threshold generation function in IntL
    g = IntL.GenThd(y,p);
    q = Q;
else
    % else generate threshold with optimization algorithm
    Gmax = IntL.Gmax;
    Gmin = IntL.Gmin;
    Qmax = IntP.Qmax;
    Qmin = IntP.Qmin;
    if Gmin==Gmax
        fun = @(x) abs(logsum(PosRat(IntL,IntP,y,[Gmax;x]),"all")-log(Nsam)-log(p));
        if isinf(Qmax) && isinf(Qmin)
            for i = 1:maxi
                [q_tep,err] = fminsearch(fun,Q,options);
                if err<err_tol
                    break
                end
            end
        elseif isinf(Qmax)
            for i = 1:maxi
                [q_tep,err] = fmincon(fun,Q,[],[],[],[],Qmin,[],[],options);
                if err<err_tol
                    break
                else
                    Q = max(q_tep,Q);
                end
            end
        elseif isinf(Qmin)
            for i = 1:maxi
                [q_tep,err] = fmincon(fun,Q,[],[],[],[],[],Qmax,[],options);
                if err<err_tol
                    break
                else
                    Q = min(q_tep,Q);
                end
            end
        else
            for i = 1:maxi
                [q_tep,err] = fminbnd(fun,Qmin,Qmax,options);
                if err<err_tol
                    break
                end
            end
        end
        q = q_tep;
        g = G;
    elseif Qmin==Qmax
        fun = @(x) abs(logmean(ResRat(IntL,IntP,y,[x;Qmax]),"all")-log(p));
        if isinf(Gmax) && isinf(Gmin)
            for i = 1:maxi
                [g_tep,err] = fminsearch(fun,G,options);
                if err<err_tol
                    break
                end
            end
        elseif isinf(Gmax)
            for i = 1:maxi
                [g_tep,err] = fmincon(fun,G,[],[],[],[],Gmin,[],[],options);
                if err<err_tol
                    break
                else
                    G = max(g_tep,G);
                end
            end
        elseif isinf(Gmin)
            for i = 1:maxi
                [g_tep,err] = fmincon(fun,G,[],[],[],[],[],Gmax,[],options);
                if err<err_tol
                    break
                else
                    G = min(g_tep,G);
                end
            end
            
        else
            for i = 1:maxi
                [g_tep,err] = fminbnd(fun,Gmin,Gmax,options);
                if err<err_tol
                    break
                end
            end
        end
        g = g_tep;
        q = Q;
    else
        fun = @(x) abs(logsum(PosRat(IntL,IntP,y,x),"all")-log(Nsam)-log(p));
        if IntL.FlgInc && IntP.FlgInc
            gqmin = [G;Q];
            gqmax = [Gmax;Qmax];
        elseif IntL.FlgInc && ~IntP.FlgInc
            gqmin = [G;Qmin];
            gqmax = [Gmax;Q];
        elseif ~IntL.FlgInc && IntP.FlgInc
            gqmin = [Gmin;Q];
            gqmax = [G;Qmax];
        else
            gqmin = [Gmin;Qmin];
            gqmax = [G;Q];
        end
        gq0 = rand(2,1).*(gqmax-gqmin)+gqmin;
        gq = fmincon(fun,gq0,[],[],[],[],gqmin,gqmax);
        g = gq(1);
        q = gq(2);
    end
end