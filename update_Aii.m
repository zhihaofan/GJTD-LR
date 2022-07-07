function [Aii] = update_Aii(par_set, B1, Lamda1, ci)
mu_1 = par_set.mu_1;
mu_2 = par_set.mu_2 * 2;
class_Lamda1 = squeeze(Lamda1(ci, :, :, :));

Aii = soft(B1-class_Lamda1/(2*mu_1),mu_1/(2*mu_2));

return;