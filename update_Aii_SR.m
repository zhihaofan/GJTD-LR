function [Aii] = update_Aii_SR(par_set, B1, Lamda1, ci, TrainLabel)
mu_1 = par_set.mu_1;
mu_2 = par_set.mu_2 * 2;
class_Lamda1 = Lamda1{ci};

Aii = soft(B1-class_Lamda1(TrainLabel == ci, :, :)/(2*mu_1),mu_1/(2*mu_2));

return;