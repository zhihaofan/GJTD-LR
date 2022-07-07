function [M] = merge_m(Dh, A_iter, ci, idx, Fish_ipts)
M1 = tprod(Dh, A_iter(:, idx, :));
M2 = tprod(Dh(:, Fish_ipts.dls == ci, :), A_iter(Fish_ipts.dls == ci, idx, :));

M = (M1 + M2) / 2;
end