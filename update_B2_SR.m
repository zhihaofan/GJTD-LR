function [B2] = update_B2_SR(par_set, Dja, Aija, EI, H1, Lamada2, ci, TrainLabel)

mu_1 = par_set.mu_1;
mu_2 = par_set.mu_2;
S = size(Dja, 3);

ifft_Dija = fft(Dja, [], 3);
ifft_EI = fft(EI, [], 3);
ifft_H1 = fft(H1, [], 3);
ifft_Aija = fft(Aija, [], 3);
ifft_B2 = zeros(size(Aija));
select_Lamada1 = Lamada2{ci}(TrainLabel ~= ci, :, :);

for i=1:S
    B2_1 = ifft_Dija(:, :, i)' * ifft_Dija(:, :, i) + mu_1 * ifft_EI(:, :, i);
    B2_2 = ifft_Dija(:, :, i)' * ifft_H1(:, :, i) + mu_1 * (ifft_Aija(:, :, i) + select_Lamada1(:, :, i));
    ifft_B2(:, :, i) = pinv(B2_1) * B2_2;
end

B2 = ifft(ifft_B2, [], 3);
return;