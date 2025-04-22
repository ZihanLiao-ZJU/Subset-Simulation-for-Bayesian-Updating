function [lograt,c] = ResRat(intBay,y)
% reshape y
y_tmp = y(1:end,:);
intBay.Y = y_tmp(3:end,:);
Li = intBay.EvlLKF;
pi = intBay.EvlPDF;
lograt_tmp = Li+pi-y_tmp(1,:)-y_tmp(4,:);
c = max(lograt_tmp,[],"all");
lograt_tmp = lograt_tmp-c;
lograt_tmp(isnan(lograt_tmp)) = -inf;
% reshape lograt
lograt = zeros(size(y,2:3));
lograt(:) = lograt_tmp;
end