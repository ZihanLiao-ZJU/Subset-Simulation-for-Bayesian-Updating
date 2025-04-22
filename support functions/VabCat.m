function varargout = VabCat(dim,varargin)
% Shrink the input variables
% Syntax:
% ----------------------------------------------------------------------
% varargout = VabSrk(Nl,varargin)
% ----------------------------------------------------------------------
% INPUTS:
% Nl        : the size shrink to
% varargin  : variables to be shrunk
% ----------------------------------------------------------------------
% OUTPUTS:
% varargout : shrunk variables
% ----------------------------------------------------------------------
Nout = nargout;
varargout = cell(Nout,1);
for iout = 1:nargout
    varargout{iout} = cat(dim,varargin{2*(iout-1)+1},varargin{2*iout});
end
end