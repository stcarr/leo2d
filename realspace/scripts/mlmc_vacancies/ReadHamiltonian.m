%
% Read Hamiltonian
%
load ml_mos2_mlmc_matrix.dat;
H = spconvert(ml_mos2_mlmc_matrix);

%%
% Compute eigenvalues
fprintf(1,'[BGN] compute eigenvalues\n');

D = eig(full(H));

fprintf(1,'[END] compute eigenvalues\n');

%%
% Estimate density of states
[mH,nH] = size(H);
de = 0.01;
regeps = 0.01;
[dos,eng] = compute_dos(D,de,regesp);
dos = dos/(mH*nH);
emin = min(D); emax = max(D);
[dos1,eng1]=ksdensity(D,eng,'kernel','normal','function','pdf','bandwidth',0.01,'support',[-1,0.8]);

figure(1);
plot(eng,dos); hold on;


figure(2);
plot(eng1,dos1); hold on;
