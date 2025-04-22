function sample=RS
minmax = [-25,15;-15,15;0,13;0,13;6,14;6,14];
sample=cell(1,6);

for i=1:6
    xtmp=[];
    N=0;
    Nsam=20000;
    xmin=minmax(i,1);
    xmax=minmax(i,2);
    g=1/(xmax-xmin);

    p = rand(1,Nsam);
    while N<=100000
        x=rand(1,Nsam).*(xmax-xmin)+xmin;
        y=MargPDF(x,i);
        %     c=max(y)/g;
        c=10;
        y=y/(c*g);
        ind_rej= y< p;
        Nacp = sum(~ind_rej);
        N=N+Nacp;
        xtmp=[xtmp x(~ind_rej)];
    end
    sample(i)={xtmp};
end
end


function pdv = MargPDF(theta,idim)
beta_a=1;alpha_a=1;c_a=10;
beta_b=1;alpha_b=1;c_b=-10;
mu_c=10;sd_c=1;
mu_d=-10;sd_d=1;
beta_i=1;alpha_i=1;c_i=10;
mu_i=10;sd_i=1;
if idim==1
    % mix log Gamma
    Log_g_a = alpha_a*(theta-c_a) - exp(theta-c_a)/beta_a - log(gamma(alpha_a)) - alpha_a*log(beta_a);
    Log_g_b = alpha_b*(theta-c_b) - exp(theta-c_b)/beta_b - log(gamma(alpha_b)) - alpha_b*log(beta_b);
    Log_g_ab = max(Log_g_a,Log_g_b);
    pdv = exp(log(1/2) + Log_g_ab + log(exp(Log_g_a-Log_g_ab) + exp(Log_g_b-Log_g_ab)));
elseif idim==2
    % mix Normal
    Log_n_c = -1/2 * ((theta-mu_c)/sd_c).^2 -log(sqrt(2*pi)) - log(sd_c);
    Log_n_d = -1/2 * ((theta-mu_d)/sd_d).^2 -log(sqrt(2*pi)) - log(sd_d);
    Log_n_cd = max(Log_n_c,Log_n_d);
    pdv = exp(log(1/2) + Log_n_cd + log(exp(Log_n_c-Log_n_cd) + exp(Log_n_d-Log_n_cd)));
elseif idim>=3 && idim<=4
    pdv = exp(alpha_a*(theta-c_i) - exp(theta-c_i)/beta_a - log(gamma(alpha_i)) - alpha_a*log(beta_i));
elseif idim>=5 && idim<=6
    pdv = exp(-1/2 * ((theta-mu_i)/sd_i).^2 -log(sqrt(2*pi)) - log(sd_i));
end
end