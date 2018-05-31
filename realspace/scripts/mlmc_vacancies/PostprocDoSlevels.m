%
% Postprocess the simulation 
% Store the averaged CCCM 
% Plot Expected value of CCCM and Variance
%
    clear all;
% 
    plot2D=0;
    col = ['r','b','g','y','m','r','b','g','y','m'];
    observable_name = 'DoS(dE) measure: ';
    picsformat = 'png';
    workdir = pwd;
    ixdir = strfind(workdir,'/');
    dirname = workdir(ixdir(end-1)+1:ixdir(end)-1);
    emfig = ['EDoS_',dirname];
    varfig= ['VarDoS_',dirname];
    stdfig= ['StdDoS_',dirname];
    
    if exist('./problem_definition.mat','file') == 2
        %
        % Config file stored during generating (newer versions of problem
        % definitions
        %
        load ./problem_definition.mat;
        nlev = max_lev;
        nslev = M_v(1:nlev);
        p = poly_order;
        fprintf(1,'==> Problem definition                       <==\n');
        fprintf(1,'   Number of levels:      %d  of %d\n',nlev,Nr_ext+1);
        fprintf(1,'   Polynomial order:      %d \n', p );
        fprintf(1,'   Defect probabilities:  Mo: %4.3f  X_A: %4.3f  X_B: %4.3f \n',pVacs)
        fprintf(1,'   Levels:\n');
        for i=1:nlev
            fprintf(1,'   l = %d size = %d x %d  # samples = %d \n',i,...
                N_v(i),N_v(i),nslev(i));
        end
        anns1 = ['p = Mo:',num2str(pVacs(1)),' X_A:',num2str(pVacs(2))];
        anns2 = ['Ns = '];
        anns3 = ['L = '];
        anns4 = ['poly = ',int2str(poly_order)];
        for ilev=1:nlev
            anns2=[anns2,int2str(nslev(ilev)),' '];
            anns3=[anns3,int2str(N_v(ilev)),'x',int2str(N_v(ilev)),' '];
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
    D = dir('./*.bin');
    nd = length(D);
    maxlevel = 0;
    for id = 1:nd
        strproc=D(id).name;
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
         if strfind(strproc,'_ENERGIES')
            idxlevelENERGY(level)=id;
        end
    end
    
    
    fprintf(1,'Number of levels: %d\n',maxlevel);
    fprintf(1,'==> Assembling the MCML estimator \n');
    
    
    %
    % Start processing data to compute MLMC estimator
    M = zeros(p,1);
    V = zeros(p,1);
   
    ilev = 1;
    while ilev <= maxlevel
        fprintf(1,'    Level: %d \n',ilev);
        
        if ilev < maxlevel
            id = idxlevelE(ilev);
            file_name = file{id}.name;
            fM = fopen(file_name,'r');       % open the E(DoS)
            Ml = fread(fM,[p 1],'double');   % p x 1 matrix load
            ML{ilev}=Ml;
            fclose(fM);
            fprintf(1,'      File: %s \n',file_name);
            
            id = idxlevelDE(ilev);
            file_name = file{id}.name;
            fDM = fopen(file_name,'r');       % open the DE(DoS)
            DMl = fread(fDM,[p 1],'double');  % p x 1 matrix load
            fclose(fDM);
            fprintf(1,'      File: %s \n',file_name);
            
            id = idxlevelV(ilev);
            file_name = file{id}.name;
            fV  = fopen(file_name,'r');        % open the Var(DoS)
            M2l= fread(fV,[p,1],'double');     % p x 1 matrix load
            VML{ilev} = M2l - Ml.^2;          
            fclose(fV);
            fprintf(1,'      File: %s \n',file_name);
            
            id = idxlevelDV(ilev);
            file_name = file{id}.name;
            fDV  = fopen(file_name,'r');       % open the DVar(M)
            DM2l= fread(fDV,[p,1],'double');   % p x 1 matrix load
            fclose(fDV);
            fprintf(1,'      File: %s \n',file_name);
 
            %
            % Accumulate estimator of Expected value of E[M]
            if ilev == 1
                M = M + Ml + DMl;
                V = V + (M2l - Ml.^2)/(M_v(ilev)-1);
            else
                M = M + DMl;
                V = V + (DM2l - DMl.^2)/(M_v(ilev)-1);
            end
        else
            id = idxlevelE(ilev);
            file_name = file{id}.name;
            fM = fopen(file_name,'r');       % open the E(DoS)
            Ml = fread(fM,[p 1],'double');   % p x 1 matrix load
            ML{ilev}=Ml;
            fclose(fM);
            fprintf(1,'      File: %s \n',file_name);
            
            id = idxlevelV(ilev);
            file_name = file{id}.name;
            fV  = fopen(file_name,'r');       % open the Var(DoS)
            M2l= fread(fV,[p,1],'double');    % p x 1matrix load
            VML{ilev} = (M2l - Ml.^2)/(M_v(ilev)-1);
            fclose(fV);
            fprintf(1,'      File: %s \n',file_name);
            
            id = idxlevelENERGY(ilev);
            file_name = file{id}.name;
            fE = fopen(file_name,'r');     % open the file with energy mesh
            E = fread(fE,[p 1],'double');  % p energy samps load
            fclose(fE);
            fprintf(1,'      File: %s \n',file_name);
        end
        ilev = ilev + 1;
    end

    %%
    % Plotting figures
    %
    %----------------------------------------------------------------------
    %
    % Plotting estimates of expected values
    %
    % Plot the final value of E[DoS]
    strtitle = [observable_name, 'MLMC estimator: E_{MLMC}[DoS]'];
    ztitle = 'E[DoS(E)]';
    
    figure(1);clf;
    %leg(1)=plot(E,M/N_v(end)^2,'k','LineWidth',2); hold on;
    leg(1)=plot(E,M,'k','LineWidth',2); hold on;
    strleg{1} = 'E_{MLMC}[DoS] estimator';
    for ilev=1:maxlevel
        %leg(ilev+1)=plot(E,ML{ilev}/N_v(ilev)^2,col(ilev)); hold on;
        leg(ilev+1)=plot(E,ML{ilev},col(ilev)); hold on;

        strleg{ilev+1} = sprintf('E_l[DoS] l = %d',ilev);
    end
    xlabel('E'); ylabel('E[DoS(E)]');
    strtitle = [observable_name,' E[DoS(E)]']; 
    title(strtitle);
    legend(leg,strleg);
    annotation('textbox',[0.15, 0.7,0.3,0.18],...
        'String',{anns1,anns2,anns3,anns4},...
        'FontSize',9,'FitBoxToText','on');
    annotation('textbox',[0.15, 0.15,0.1,0.05],...
        'String',[dirname],...
        'FontSize',9,'LineStyle','none','FitBoxToText','on');
    saveas(gcf,emfig,picsformat);
   
    figure(2);clf;
       plot(M/N_v(end)^2,E,'k','LineWidth',2); hold on;
       ylabel('E'); xlabel('E[DoS(E)]');
       strtitle = [observable_name,' E[DoS(E)] Poly order=',int2str(p)]; 
       title(strtitle);
       
 %----------------------------------------------------------------------

 %----------------------------------------------------------------------
 %
 % Plotting Variance estimates
 %
 % Plot the final value of E[M]
    strtitle = [observable_name, 'MLMC estimator: Var_{MLMC}[DoS]'];
    ztitle = 'Var[DoS(E)]';
    
%%    
    figure(3);clf;
    %leg(1)=plot(E,V,'k','LineWidth',2); hold on;
    leg(1)=semilogy(E,V,'k','LineWidth',2); hold on;

    strleg{1} = 'Var_{MLMC}[DoS] estimator';
    for ilev=1:maxlevel
        %leg(ilev+1)=plot(E,VML{ilev},col(ilev)); hold on;
        leg(ilev+1)=semilogy(E,VML{ilev},col(ilev)); hold on;
        strleg{ilev+1} = sprintf('Var_l[DoS] l = %d',ilev);
    end
    xlabel('E'); ylabel('Var[DoS(E)]');
    strtitle = [observable_name,' Var[DoS(E)]']; 
    title(strtitle);
    legend(leg,strleg);
    annotation('textbox',[0.15, 0.7,0.3,0.18],...
        'String',{anns1,anns2,anns3,anns4},...
        'FontSize',9,'FitBoxToText','on');
    annotation('textbox',[0.15, 0.15,0.1,0.05],...
        'String',[dirname],...
        'FontSize',9,'LineStyle','none','FitBoxToText','on');
     saveas(gcf,varfig,picsformat);
     
    figure(4);clf;
    leg(1)=semilogy(E,sqrt(V),'k','LineWidth',2); hold on;
    strleg{1} = 'std_{MLMC}[DoS] estimator';
    for ilev=1:maxlevel
        leg(ilev+1)=semilogy(E,sqrt(VML{ilev}),col(ilev)); hold on;
        strleg{ilev+1} = sprintf('std_l[M] l = %d',ilev);
    end
    xlabel('E'); ylabel('std[DoS(E)]');
    strtitle = [observable_name,' std[DoS(E)]']; 
    title(strtitle);
    legend(leg,strleg);
    annotation('textbox',[0.15, 0.7,0.3,0.18],...
        'String',{anns1,anns2,anns3,anns4},...
        'FontSize',9,'FitBoxToText','on');
    annotation('textbox',[0.15, 0.15,0.1,0.05],...
        'String',[dirname],...
        'FontSize',9,'LineStyle','none','FitBoxToText','on');
     saveas(gcf,stdfig,picsformat);   
%     figure(8);clf;
%        plot(varccmdiag,E,'k','LineWidth',2); hold on;
%        ylabel('E'); xlabel('Var[M(E,E)]');
%        strtitle = [observable_name,' Var[M(E,E)] Poly order=',int2str(p)]; 
%        title(strtitle);
       
 %----------------------------------------------------------------------      
    