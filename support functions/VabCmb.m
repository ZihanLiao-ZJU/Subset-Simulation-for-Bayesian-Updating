function varargout = VabCmb(varargin)
Nout = nargout;
varargout = cell(Nout,1);
for iout = 1:nargout
    if ~iscell(varargin{iout})
        Ndim = size(varargin{iout},1);
        varargout{iout} = varargin{iout}(1:Ndim,:);
    else
        Ncel = length(varargin{iout});
        Ndim = zeros(Ncel,1);
        var_tep = varargin{iout};
        for icel = 1:Ncel
            Ndim(icel) = size(var_tep{icel},1);
        end
        varargout{iout} = [];
        ind_ept = Ndim==0;
        Ndim = Ndim(~ind_ept);
        if all(Ndim(1)==Ndim)
            Ndim = Ndim(1);
            for icel = 1:Ncel
                if ~ind_ept(icel)
                    varargout{iout} = cat(2,varargout{iout},var_tep{icel}(1:Ndim,:));
                end
            end
        else
            for icel = 1:Ncel
                if ~ind_ept(icel)
                    varargout{iout} = cat(2,varargout{iout},var_tep{icel}(:)');
                end
            end
        end
    end
end
end