%
% Postprocess the simulation 
% Store the averaged CCCM 
% Plot Expected value of CCCM and Variance
%
    clear all;
% 
    p = 500;                                        % poly_order
    
    D = dir('./');
    nd = length(D);
    maxlevel = 0;
    for id = 1:nd
        strproc=D{id}.name;
        is = strfind(strproc,'L');
        level = str2num(strproc(is+1));
        maxlevel = max(maxlevel,level);
        file{id}.name = strproc;
        file{id}.level = level;
        if strfind(strproc,'_E_')
            idxlevelE(level)=id;
        end
        if strfind(strproc,'_DE_')
            idxlevelDE(level)=id;
        end
        if strfind(strproc,'_V_')
            idxlevelV(level)=id;
        end
        if strfind(strproc,'_DV_')
            idxlevelDV(level)=id;
        end
    end
    fprintf(1,'Number of levels: %d\n',maxlevel);
    fprintf(1,'==> Assembling the MCML estimator \n');
    M = zeros(p,p);
    V = zeros(p,p);
    for l=1:maxlevel
        idx = idxlevelE(l);
        fname = file(idx).name;
        fM = fopen(fname,'r');   % open the E(M)_XX file
        MatIn = fread(fM,[p p],'double');
        if l == 1
            M = M + MatIn;
        else
            
        end
        
        fV = fopen('ml_mos2_mlmc_L3_V_M_xx.bin','r');   % open the Var(M)
    end
    fE = fopen('ml_mos2_mlmc_L3_ENERGIES.bin','r'); % open the Energy file
    M = fread(fM,[p p],'double');                   % p x p matrix load
    M2= fread(fV,[p,p],'double');
    E = fread(fE,[p 1],'double');                   % p energy samps load
    fclose(fM);
    fclose(fV);
    fclose(fE);

    figure(1);clf;
    surf(E,E,M,'EdgeColor','none'); hold on; % simple plot
    view(2);
    axis square;
    xlabel('Energy E_1'); 
    ylabel('Energy E_2');
    zlabel('E[M(E_1,E_2)]');
    title('Current-current correlation measure');
   
    VM = M2 - M.^2;
    figure(2);clf;
    surf(E,E,VM,'EdgeColor','none'); hold on; % simple plot
    view(2);
    axis square;
    xlabel('Energy E_1'); 
    ylabel('Energy E_2');
    zlabel('Var[M(E_1,E_2)]');
    title('Variance of the Current-current correlation measure');    
    
    ccmdiag = diag(M);
    ccmdiag2 = diag(VM);
    varccmdiag = ccmdiag2-ccmdiag.^2;
    figure(3);clf;
    subplot(2,1,1);
        plot(E,ccmdiag); hold on;
        xlabel('E');
        ylabel('E[M(E,E)]');
        title('Current-current correlation measure M(E,E)');
    subplot(2,1,2);
        plot(E,varccmdiag); hold on;
        xlabel('E');
        ylabel('Var[M(E,E)]');
        title('Variance Current-current correlation measure M(E,E)');