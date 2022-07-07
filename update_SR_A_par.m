function [A, Init_par] = update_SR_A_par(Init_par, X, D, A_iter, Fish_ipts, seg, s_idx, A_l, par)
par_set.mu_1 = 0.0001;
par_set.mu_2 = 0.5;
k3 = size(X, 3);

A = A_iter;
iii = 0;
while iii < 3
    for ci=1:length(seg)-1
        idx = s_idx(seg(ci)+1:seg(ci+1));
        Zi=X(:, idx, :);
        Di=D(:,Fish_ipts.dls ==ci,:);
        EI = zeros(size(Di, 2), size(Di, 2), size(Di, 3));
        EI(:, :, 1) = eye(size(Di, 2));
        Aija = A(Fish_ipts.dls ~= ci, idx,:);
        Dja = D(:, Fish_ipts.dls ~=ci,:);
        H = Zi - tprod(Dja, Aija);
        Aii = A(Fish_ipts.dls == ci, idx,:);
        B1 = update_B1_SR(par_set, D, EI, H, Zi, Aii, Init_par.Lamada1, ci, Fish_ipts.dls, idx);
        Aii = update_Aii_SR(par_set, B1, Init_par.Lamada1, ci, Fish_ipts.dls);
        A(Fish_ipts.dls == ci, idx, :) = Aii;
        Init_par.Lamada1{ci}(Fish_ipts.dls == ci, :, :) = Init_par.Lamada1{ci}(Fish_ipts.dls == ci, :, :) + par_set.mu_1 * (Aii - B1);
        EI1 = zeros(size(Dja, 2), size(Dja, 2), size(Dja, 3));
        EI1(:, :, 1) = eye(size(Dja, 2));
        H1 = Zi - tprod(Dja,A(Fish_ipts.dls ~= ci, idx, :));
        B2 = update_B2_SR(par_set, Dja, Aija, EI1, H1, Init_par.Lamada1, ci, Fish_ipts.dls);
        Aija = update_Aij_SR(par_set, B2, Init_par.Lamada1, ci, Fish_ipts.dls);
        A(Fish_ipts.dls ~= ci, idx, :) = Aija;
        Init_par.Lamada1{ci}(Fish_ipts.dls ~= ci, :, :) = Init_par.Lamada1{ci}(Fish_ipts.dls ~= ci, :, :) + par_set.mu_1 * (Aija - B2);
    end
    A = (A_l + A)/2;
    iii = iii + 1;
end


return;