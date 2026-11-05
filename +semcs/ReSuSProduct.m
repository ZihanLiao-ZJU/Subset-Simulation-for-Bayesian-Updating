function out=ReSuSProduct(raw,sampler)
% Native Re-SuS product, with the final TRUE failure factor (not hybrid mass).
K=raw.Nite;cp=zeros(1,K);c=zeros(1,K);thresholds=zeros(K,2);samples=cell(1,K);values=samples;history={};
for j=1:K
 state=raw.g{j};thresholds(j,:)=[state.g_threshold,state.logpdf_threshold];
 samples{j}=reshape(raw.x{j},size(raw.x{j},1),[]);values{j}=reshape(raw.y{j}(4,:,:),1,[]);
 if j<K
  next=raw.g{j+1};qrow=3;if isfield(state,'schema'),qrow=state.schema.reference_row;end
  lp=reshape(raw.y{j}(qrow,:,:),1,[]);
  cp(j)=mean(values{j}<=next.g_threshold | lp<=next.logpdf_threshold);c(j+1)=c(j)+log(cp(j));
 end
 if j>1,history{end+1}=struct('initial',raw.g{j},'final',raw.sampled_final_G{j});end
end
cp(K)=mean(values{K}<=0);logpf=c(K)+log(cp(K));
out=struct('estimator','product','conditional_probabilities',cp,'sampler',sampler, ...
 'profile',raw.g{1}.profile,'sampler_history',{history},'samples',{samples},'values',{values}, ...
 'thresholds',thresholds,'log_normalizers',c,'sample_counts',raw.Nsam(:)', ...
 'n_calls',sum(raw.Ncal),'converged',raw.flg_cvg,'log_pf',logpf,'pf',exp(logpf),'raw',raw);
end
