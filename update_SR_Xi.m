function [Yi] = update_SR_Xi(Mi)
eta = 1e-3;
% eta = par;

mu = 1e-3;
Yi = Log_prox_tnn( Mi, eta/2/mu);
% Yi = Log_prox_tnn( Mi, par);
% Yi = Mi;
return;