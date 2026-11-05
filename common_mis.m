function [value,parts]=common_mis(raw,C)
% Optional cached balance MIS. Query the actual model's declared support.
% Only sampled proposal levels are included; no target evaluation or RNG.
if nargin<2,post=NumSim(raw);C=post.C;end
K=numel(raw.y);C=C(:);
assert(numel(C)==K,'SuSCommonMIS:Normalizer','One normalizer per sampled proposal is required.');
target=cell(K,1);dens=cell(K,1);model=raw.intBay;
for k=1:K
    yi=raw.y{k};target{k}=yi(3,:)+yi(4,:);dens{k}=zeros(K,numel(target{k}));
    for j=1:K
        model.G=raw.g{j};model.Y=yi(3:end,:);
        dens{k}(j,:)=model.EvlLKF+model.EvlPDF-C(j);
    end
end
[value,parts]=semcs.balance_log_integral(target,dens);
end
