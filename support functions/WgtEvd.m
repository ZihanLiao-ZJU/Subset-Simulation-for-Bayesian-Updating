function w_evd = WgtEvd(g,y)
% initialization
Nite = length(y);
w_evd = cell(Nite,1);
Nsam = zeros(Nite,1);
Nsze = zeros(Nite,2);
for ite = 1:Nite
    Nsze(ite,:) = size(y{ite},2:3);
    Nsam(ite) = prod(Nsze(ite,:));
    w_evd{ite} = zeros(Nsze(ite,:));
end

for ite = 1:Nite-1
    Li = g{ite}(3);
    Li1 = g{ite+1}(3);
    w_evd_tmp = min(logminus(y{ite}(3,:),Li),logminus(Li1,Li));
    w_evd{ite}(:) = w_evd_tmp(:);
end
Li = g{Nite}(3);
w_evd_tmp = logminus(y{Nite}(3,:),Li);
w_evd{Nite}(:) = w_evd_tmp(:);
end