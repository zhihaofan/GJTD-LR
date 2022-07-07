function [B1] = update_B1(par_set, D, EI, H, Zi, Aii, Lamada1, ci, TrainLabel)

mu_1 = par_set.mu_1;
mu_2 = par_set.mu_2;
mu = 0.001;
Di = D(:,TrainLabel == ci,:);
S = size(D, 3);
Xi = Zi;

ifft_H = fft(H, [], 3);
ifft_Di = fft(Di,[],3);
ifft_EI = fft(EI, [], 3);
ifft_Xi = fft(Xi,[],3);
ifft_Aii = fft(Aii,[],3);
ifft_B1 = zeros(size(Aii));
ifft_Lamada1 = fft(Lamada1, [], 3);
select_Lamada1 = squeeze(Lamada1(ci, :, :, :));
for i=1:S
    iB1_1 = ifft_Di(:, :, i)' * ifft_Di(:, :, i) + mu_1 * ifft_EI(:, :, i);
    iB1_2 = ifft_Di(:, :, i)' * ifft_H(:, :, i) + ifft_Di(:, :, i)' * ifft_Xi(:, :, i) + ...
        mu_1 * (ifft_Aii(:, :, i) + select_Lamada1(:, :, i)/(2*mu_1));
    ifft_B1(:, :, i) =  pinv(iB1_1) * iB1_2;
end

B1 = ifft(ifft_B1, [], 3);

return