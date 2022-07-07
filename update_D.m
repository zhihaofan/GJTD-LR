function [new_D] = update_D(par_set, D, Zi, A, ci, TrainLabel, Fish_ipts)
mu_1 = par_set.mu_1 * 100;
mu_2 = par_set.mu_2;
Xi = Zi;
Aii = A(TrainLabel == ci,Fish_ipts.trls == ci,:);
Aija = A(TrainLabel ~= ci,Fish_ipts.trls == ci,:);
Di = D(:,TrainLabel ==ci,:);
Dja = D(:,TrainLabel ~=ci,:);
new_D = zeros(size(D));

ifft_Xi = fft(Xi,[],3);
ifft_Aii = fft(Aii,[],3);
ifft_Aija = fft(Aija,[],3);
ifft_Di = fft(Di,[],3);
ifft_Dja = fft(Dja,[],3);

S = size(Zi, 3);
for i=1:S
    iDi_1 = ifft_Xi(:, :, i) * ifft_Aii(:, :, i)' + ifft_Dja(:, :, i) * ifft_Aija(:, :, i) * ifft_Aii(:, :, i)';
    iDi_2 = ifft_Aii(:, :, i) * ifft_Aii(:, :, i)';
    ifft_Di(:, :, i) = mu_1 * iDi_1 * pinv(iDi_2);
    
end
new_Di = ifft(ifft_Di, [], 3);
new_D(:,TrainLabel ~=ci,:) = Dja;
new_D(:,TrainLabel ==ci,:) = new_Di;

return;