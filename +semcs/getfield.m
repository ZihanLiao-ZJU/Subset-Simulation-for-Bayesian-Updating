function value=getfield(g,path,default)
parts=strsplit(path,'.');value=g;
for j=1:numel(parts)
 if ~isstruct(value)||~isfield(value,parts{j}),value=default;return;end
 value=value.(parts{j});
end
end
