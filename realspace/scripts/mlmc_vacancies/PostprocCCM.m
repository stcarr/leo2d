%
% Postprocess the simulation 
% Store the averaged CCCM 
% Plot Expected value of CCCM and Variance
%
    clear all;
% 
    if exist('../problem_definition.mat','file') == 2
        %
        % Config file stored during generating (newer versions of problem
        % definitions
        %
        load problem_definition.mat;
        nlev = max_lev;
        nslev = M_v(1:nlev);
        p = poly_order;
        fprintf(1,'==> Problem definition                       <==\n');
        fprintf(1,'   Number of levels:      %d  of %d\n',nlev,Nr_ext);
        fprintf(1,'   Polynomial order:      %d \n', p );
        fprintf(1,'   Defect probabilities:  Mo: %4.3f  X_A: %4.3f  X_B: %4.3f \n',pVacs)
        fprintf(1,'   Levels:\n');
        for i=1:nlev
            fprintf(1,'   l = %d size = %d x %d  # samples = %d \n',i,...
                N_v(i),N_v(i),nslev(i));
        end
    else
        % 
        % Config file not stored (older versions)
        % Set the values manually
        fprintf(1,'WARNING: problem_definition.mat not found \n');
        fprintf(1,'         check that the set values are correct\n');
        %
        nlev = 3;                                       % number of levels
        nslev=[100,20,3];                               % samples per level
        p = 500;                                        % poly_order
        fprintf(1,'==> Problem definition                       <==\n');
        fprintf(1,'   Number of levels:  %d  of %d\n',nlev,nlev);
        fprintf(1,'   Polynomial order:  %d \n', p );
        fprintf(1,'   Levels:\n');
        for i=1:nlev
            fprintf(1,'   l = %d size = unknwon   # samples = %d \n',i,...
                nslev(i));
        end
    end
    %
    % Start processing data to compute MLMC estimator
    M = zeros(p,p);
    V = zeros(p,p);
   
    
    for i = 1:nlev

        fM = fopen(['ml_mos2_mlmc_L' num2str(i) '_E_M_xx.bin'],'r');   % open the E(M)
        Ml = fread(fM,[p p],'double');                  % p x p matrix load
        ML{i}=Ml;
        fclose(fM);
        if (i < max_lev)
            fDM = fopen(['ml_mos2_mlmc_L' num2str(i) '_DE_M_xx.bin'],'r');   % open the E(M)
            DMl = fread(fDM,[p p],'double');                  % p x p matrix load
            fclose(fDM);
        end

        %fV = fopen('ml_mos2_mlmc_L1_V_M_xx.bin','r');  % open the Var(M)
        %M2l= fread(fV,[p,p],'double');
        %close(fV);

        if (i == 1)
            M = M + Ml + DMl;
        elseif (i < max_lev)
            M = M + DMl;
        end
    end
    
    fE = fopen('ml_mos2_mlmc_L3_ENERGIES.bin','r'); % open the Energy file
    E = fread(fE,[p 1],'double');                   % p energy samps load
    fclose(fE);

    %
    % Plotting figures
    %
    figure(1);clf;
    surf(E,E,M,'EdgeColor','none'); hold on; % simple plot
    view(2);
    axis square;
    xlabel('Energy E_1'); 
    ylabel('Energy E_2');
    zlabel('E[M(E_1,E_2)]');
    title('Current-current correlation measure');
   
    figure(2);clf;
    subplot(2,2,1);
       surf(E,E,ML{1},'EdgeColor','none'); hold on; % simple plot
       view(2);
       axis square;
       xlabel('Energy E_1'); 
       ylabel('Energy E_2');
       zlabel('E[M(E_1,E_2)]');
       title('Current-current correlation measure Level 1');  
    subplot(2,2,2);   
    surf(E,E,ML{2},'EdgeColor','none'); hold on; % simple plot
       view(2);
       axis square;
       xlabel('Energy E_1'); 
       ylabel('Energy E_2');
       zlabel('E[M(E_1,E_2)]');
       title('Current-current correlation measure Level 2'); 
    subplot(2,2,3);   
    surf(E,E,ML{3},'EdgeColor','none'); hold on; % simple plot
       view(2);
       axis square;
       xlabel('Energy E_1'); 
       ylabel('Energy E_2');
       zlabel('E[M(E_1,E_2)]');
       title('Current-current correlation measure Level 3'); 
    subplot(2,2,4);   
    surf(E,E,M,'EdgeColor','none'); hold on; % simple plot
       view(2);
       axis square;
       xlabel('Energy E_1'); 
       ylabel('Energy E_2');
       zlabel('E[M(E_1,E_2)]');
       title('Current-current correlation measure MLMC estimator');   
       
    
    ccmdiag = diag(M);
    for i=1:nlev
        Ml = ML{i};
        CCMdiag{i} = diag(Ml);
    end
    figure(3);clf;
    
        leg(1)=plot(E,ccmdiag,'k','LineWidth',2); hold on;
        
        for i = 1:max_lev
            leg(1+i)=plot(E,CCMdiag{i}); hold on;
        end

        xlabel('E');
        ylabel('E[M(E,E)]');
        title('Current-current correlation measure M(E,E)');
        legend(leg,'MLMC','Lev 1','Lev 2','Lev 3','Lev 4', 'Lev 5')
        
    figure(4);clf;
       plot(ccmdiag,E,'k','LineWidth',2); hold on;
       ylabel('E');
       xlabel('E[M(E,E)]');
       title('Current-current correlation measure M(E,E)');
    