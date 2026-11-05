function g=setfield(g,path,value)
parts=strsplit(path,'.');g=put(g,parts,value);
end
function g=put(g,p,value)
if ~isstruct(g),g=struct();end
if numel(p)==1,g.(p{1})=value;return;end
if ~isfield(g,p{1}),g.(p{1})=struct();end
g.(p{1})=put(g.(p{1}),p(2:end),value);
end
