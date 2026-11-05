function [value,contributions]=balance_log_integral(logtargets,crossdens)
% Balance MIS of cached targets/densities; zero target evaluations and RNG.
K=numel(logtargets);N=cellfun(@numel,logtargets);contributions=zeros(K,1);
for k=1:K
 density=crossdens{k};assert(isequal(size(density),[K,N(k)]),'semcs:MISShape','Cross density must be K by actual sample count.');
 denominator=logsum(density+log(N(:)),1);
 contributions(k)=log(N(k))+logsum(reshape(logtargets{k},1,[])-denominator,'all')-log(N(k));
end
value=logsum(contributions,'all');
end
function value=logsum(x,dim)
top=max(x,[],dim);value=top+log(sum(exp(x-top),dim));value(isinf(top))=top(isinf(top));
end
