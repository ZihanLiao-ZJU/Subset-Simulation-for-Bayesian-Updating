function name=kernel_name(s)
if isa(s,'semcs.aCS'),name='aCS';
elseif isa(s,'semcs.ES'),name='ES';
elseif isa(s,'semcs.HMC'),name='HMC';
elseif isa(s,'semcs.CW_aMH'),name='CW_aMH';
else,name=matlab.lang.makeValidName(class(s));end
end
