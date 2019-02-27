clear all;

load problem_definition.mat
p = poly_order;

tar_lev = 2;

fM = fopen(['ml_mos2_mlmc_L' num2str(tar_lev) '_E_M_xx.bin'],'r');   % open the E(M)
M = fread(fM,[p p],'double');                  % p x p matrix load
fclose(fM);

fE = fopen(['ml_mos2_mlmc_L' num2str(tar_lev) '_ENERGIES.bin'],'r');   % open the E(M)
E = fread(fE,[1 p],'double');                  % p x p matrix load
fclose(fE);

E_f = -1.5;
T = 100;
w_list = [0:.01:6];


sigma_xx = zeros(size(w_list));

V = (10^-8); % volume of sample, just a prefactor
kb = 8.6173303e-5; % Boltzman's Constant, eV / K
hbar = 6.5821e-16; % only effects 
eta_hbar = 0.03; % hbar*eta


fermi = @(x,y,z) 1./(exp((x - z)./(kb*y)) + 1);
max_E = kb*T*20;
dE = mean(diff(-E));
E_search = ceil(max_E/dE);
[~,E_f_idx] = min(abs(E-E_f));


for w_idx = 1:length(w_list)
    w = w_list(w_idx);
    w_idx/length(w_list)
    for E1_idx = max(1,E_f_idx-E_search):min(length(E),E_f_idx+E_search)
        E_n = E(E1_idx);
        
        [~,E_w_idx] = min(abs(E - (E_n + w)));
        
        for E2_idx =max(1,E_w_idx-E_search):min(length(E),E_w_idx+E_search)
            E_m = E(E2_idx);
            
            if E1_idx ~= E2_idx
                fermi_fac = fermi(E_n,T,E_f) - fermi(E_m,T,E_f);
                denom = (E_m - E_n)*((w + 1i*eta_hbar) + (E_n - E_m));

                sigma_xx(w_idx) = sigma_xx(w_idx) + 1i*(hbar/V)*M(E1_idx,E2_idx)*fermi_fac/denom;
            end
        
        end
    end
end
%%
clf
subplot(2,1,1)
hold on
plot(w_list,imag(sigma_xx),'r');
plot(w_list,real(sigma_xx),'k');
plot([0 w_list(end)],[0 0],'--k');

subplot(2,1,2)

surf(E,E,M,'EdgeColor','none')

E_min = E_f - 0;
E_max = E_f + w_list(end) + 0;
axis([E_min E_max E_min E_max 0 inf])
view(2)
colormap hot
caxis([0 10])
