function p = setDefaultParameters()
% Setting parameters
phi=0.125;
gca=4.0;
gl=2.0;
gk=8.0;
v1=-2.8;
v2=26.0;
v3=12.0;
v4=17.4;
vl=-60.0;
vca=120.0;
vk=-91.893652410657;
C=20.0;
Kappa_s=5.0;
v_s=0.0;
I_I=-60.0;
Kappa_I=0.0;
Tau_I=10.0;
p.ind_I=13;
p.ind_Kappa=14;
p.ind_Tau=17;
p.parameters_I=[phi, gca, gl, gk, v1, v2, v3, v4, vl, vca, vk, C, I_I, Kappa_I, v_s, Kappa_s, Tau_I];
p.initialguess_I=[-100,0];

end
