function [Aij] = update_Aij_SR(par_set, B2, Lamada2, ci, TrainLabel)
mu_1 = par_set.mu_1;
mu_2 = par_set.mu_2 * 2;
class_Lamada2 =Lamada2{ci};

Aij = soft(B2-class_Lamada2(TrainLabel ~= ci, :, :)/(2*mu_1),mu_1/(2*mu_2));

return;