% Set_up_samples_MoS2.m
% 
% This script samples random vacancies on an MLMC hierarchy of MoS2
% single-layer tight-binding models to be used with Stephen Carr's "LEO2D"
%
% Modification of the original set_up_samples_MoS2.m
%
% Modified: pp 06/12/17
%

%----------------------------------------------------
% FOR READING OUTPUT FROM LEO2D:
%{

    % for conductivity
    % Do something like this:
    p = 500;                                        % poly_order
    fM = fopen('ml_mos2_mlmc_L3_E_M_xx.bin','r');   % open the M_XX file
    fE = fopen('ml_mos2_mlmc_L3_ENERGIES.bin','r'); % open the Energy file
    M = fread(fM,[p p],'double');                   % p x p matrix load
    E = fread(fE,[p 1],'double');                   % p energy samps load
    fclose(fM);
    fclose(fE);

    surf(E,E,M,'EdgeColor','none'), view(2), axis square; % simple plot

    % for DoS
    % Do something like this:
    p = 500;                                        % poly_order
    fdos = fopen('ml_mos2_mlmc_L3_E_dos.bin','r');   % open the M_XX file
    fE = fopen('ml_mos2_mlmc_L3_ENERGIES.bin','r'); % open the Energy file
    dos = fread(fdos,[p 1],'double');                   % p x p matrix load
    E = fread(fE,[p 1],'double');                   % p energy samps load
    fclose(fdos);
    fclose(fE);

    hold on
    plot(E,dos), view(2), axis square; % simple plot
    

%}
%----------------------------------------------------

fprintf(1,'\n\n***\n');
fprintf(1,'*** Set up MLMC hierarchy for the tight binding model of MoS2\n');
fprintf(1,'***\n\n');
fprintf(1,'Model:\n');

% Tight-binding model parameters
nrat = 3;                           % Nr of atoms per unit cell
orbindcs = [1,1,1,1,1,2,2,2,3,3,3]; % Atom numbers associated with orbitals
fprintf(1,'  # of atoms per unit cell: %d \n',nrat);
fprintf(1,'  atom numbers\n');
fprintf(1,'      [');
for i=1:length(orbindcs)
    fprintf(1,' %d ',orbindcs(i));
end
fprintf(1,']\n');

% Probability of vacancy in atom sites of different types:
pVacs = [0.0;0.05;0.05]; % M, X_A, X_B

% Parameters controlling the hierarchy
Nbase  = 2; % Side of smallest super cell
Nr_ext = 2;  % Number of extensions (by doubling) from base to largest cell
% The following vector prescribes number of samples, M, per level.
M_v     = [100,20,3]; % Length should be Nr_ext+1
% The choice M=100,20,3 is small for the smaller super cell sizes.  
% Number of samples on different super cell sizes must in practice be
% adjusted to parameters estimated from runs.

% Polynomial order for the Chebyshev fitting to the current-current measure
poly_order = 80;

% K sampling
k_sampling = 1;
k_grid = 1;

% Default out and scratch/temp directories
outdir = '../out';
tempdir = '../temp';

% Default solver settings
dos_on = 0;
cond_on = 1;

choice=-1;
while choice~=0
    fprintf(1,'==> Menu <==\n');
    fprintf(1,' 1 ... change vacancy probabilities \n');
    fprintf(1,'         pVacs = Mo: %4.3f  X_A: %4.3f   X_B: %4.3f\n',pVacs);
    fprintf(1,' 2 ... change the side of the smallest super cell\n');
    fprintf(1,'         Nbase = %d \n',Nbase);
    fprintf(1,' 3 ... change number of extensions by doubling the base cell\n');
    fprintf(1,'         Nr_ext = %d \n',Nr_ext);
    fprintf(1,' 4 ... change number of samples per level \n');
    for i = 1:Nr_ext+1
        Ns0=Nbase*2^(i-1);
        fprintf(1,'         level: %d  size: %d x %d   # samples: %5d \n', ...
            i,Ns0,Ns0,M_v(i));
    end
    fprintf(1,' 5 ... change polynomial order \n');
    fprintf(1,'         poly_order = %d \n',poly_order);
    fprintf(1,' 6 ... change k-point sampling (1/0 = on/off) \n');
    fprintf(1,'         k_sampling = %d \n',k_sampling);    
    fprintf(1,' 7 ... change k-grid at highest level \n');
    fprintf(1,'         k_grid = %d \n',k_grid);
    fprintf(1,' 8 ... change output and scratch directories \n');
    fprintf(1,['         outdir = ',outdir,'\n']);
    fprintf(1,['         tempdir = ',tempdir,'\n']);
    fprintf(1,' 9 ... change DoS and conductivity settings \n');
    fprintf(1,'         DoS = %d \n',dos_on);
    fprintf(1,'         Cond = %d \n',cond_on);    
    fprintf(1,' 0 ... proceed to generating config files\n');
    
    choice=input('Choice: ');
    switch choice
        case 1
            pVacs(1) = input('Probability of Mo  vacancies: ');
            pVacs(2) = input('Probability of X_A vacancies: ');
            pVacs(3) = input('Probability of X_B vacancies: ');
        case 2
            Nbase = input('Size of the smallest supercell: ');
        case 3
            Nr_ext = input('Number of extensions: ');
            while (length(M_v) < Nr_ext+1)
               M_v(end+1) = 1; 
            end
        case 4
            for i = 1:Nr_ext+1
                Ns0=Nbase*2^(i-1);
                fprintf(1,'level: %d  size: %d x %d   # samples: %5d \n', ...
                    i,Ns0,Ns0,M_v(i));                
                nsamples=input('Number of samples for this level (0=no change): ');
                if nsamples>0
                    M_v(i) = nsamples;
                end                
            end
        case 5
          poly_order = input('Polynomial order: ');
        case 6
          k_sampling = input('k sampling 1/0 (on/off): ');
        case 7
          k_grid = input('k Grid size on highest level: ');
        case 8
          outdir = input('Output directory: ','s');
          tempdir = input('Scratch directory: ','s');
        case 9
          dos_on = input('DoS (0/1 = off/on): ');
          cond_on = input('Conductivity (0/1 = off/on): ');          
    end
    
end

fprintf(1,'[BGN] Generate files with vacancies\n');
% ----------------------------------------------------------------------- %
% Generate a hierarchy of samples
% All MLMC levels except smallest super cell have a control variate, i.e. 
% are  "bi-level". 
biLev_v = [false,true(1,Nr_ext)]; % Never change this
% The control variates are based on subdividing larger super cells into
% four parts. This vector holds the side of super cell per level
N_v     = Nbase*2.^(0:Nr_ext);    % Never change this

for k=1:Nr_ext+1 
  dirstr = ['lev_',int2str(k)];
  fprintf(1,'   level: %d  size: %d x %d   # samples: %d FILE: Vac_%s\n', ...
                    k,N_v(k),N_v(k),M_v(k),dirstr);
                
  generate_samples(pVacs,biLev_v(k),N_v(k),N_v(k),nrat,...
    orbindcs,M_v(k),dirstr)
end

fprintf(1,'[END]\n\n')

answer = input('Merge the files (0/1): ');
if answer~=0
    fprintf(1,'[BGN] Merging the files');
    strprompt = sprintf('Max level to be merged (value in [1,%d]): ',Nr_ext+1);
    max_lev = input(strprompt);
    problem_dir=merge_v_files(max_lev,N_v,M_v,poly_order,k_sampling,k_grid,outdir,tempdir,dos_on,cond_on);
    fprintf(1,'Problem definition in: %s\n',problem_dir);
    fprintf(1,'[END]\n\n');
end

save ([problem_dir,'/problem_definition.mat'],'max_lev','Nbase','Nr_ext','M_v','N_v','pVacs','poly_order');


    