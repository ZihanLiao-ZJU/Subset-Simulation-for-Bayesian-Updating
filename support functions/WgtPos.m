function w_pos = WgtPos(g,y)
Nite = length(y);
w_pos = cell(Nite,1);
for ite = 1:Nite
    Nsze = size(y{ite});
    Nsze = Nsze(2:end);
    if ite<Nite
        w_pos_tep = y{ite}(3,:) + y{ite}(4,:) - y{ite}(1,:) - y{ite}(2,:) + log(double(y{ite}(3,:)<g{ite+1}(3)));
    else
        w_pos_tep = y{ite}(3,:) + y{ite}(4,:) - y{ite}(1,:) - y{ite}(2,:);
    end
    w_pos{ite} = zeros([Nsze,1]);
    w_pos{ite}(:) = w_pos_tep(:);
    w_pos{ite}(isnan(w_pos{ite})) = -inf;
end
end