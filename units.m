function [mu_e, D_e, n_0, l_0, t_0, E_0] = units(N)
% From Ute-Montijn paper: Raether-Meek Criterion
mu_e = 380*(1/N);
D_e = 1800*(1/N);
n_0 = 4.8e+14*((N/1)^2);
l_0 = 2.3e-6*(1/N);
t_0 = 3e-12*(1/N);
E_0 = 200*(N/1);
end