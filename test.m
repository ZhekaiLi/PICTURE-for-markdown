clear;
close all;
clc;


%% find tp
Ks = 30;
Phi_f = 61;
n = 0.44;
theta_i = 0.12;

ntp = (n - theta_i) * Phi_f;
dtheta = (n - theta_i);

ri = 52; % rainfall rate

syms F;
fp = Ks + Ks * ntp / F;

eq = fp - ri;
Fp = solve(eq, F);
tp = Fp/ri;

%% solve fp vs t

FList = zeros(100, 1);
fList = zeros(100, 1);
LList = zeros(100, 1);
RList = zeros(100, 1);
tList = 0:0.01:0.99;

for i = 1:100
    t = tList(i);
    if i <= tp * 100
        FList(i, 1) = ri * t;
        fList(i, 1) = ri;
        LList(i, 1) = FList(i, 1)  / dtheta;
    else
        syms F;
        eq = t == tp + 1/Ks * ((F - Fp) - Phi_f * dtheta * log((Phi_f*dtheta + F)/(Phi_f*dtheta + Fp)));
        Fi = solve(eq, F);
        
        FList(i, 1) = real(Fi(1));
        fList(i, 1) = Ks + Ks * Phi_f * dtheta / FList(i, 1);
        LList(i, 1) = FList(i, 1) / dtheta;
        RList(i, 1) = ri * t - FList(i, 1);
    end
    
end
