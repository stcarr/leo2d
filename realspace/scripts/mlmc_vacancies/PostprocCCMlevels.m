%
% Postprocess the simulation 
% Store the averaged CCCM 
% Plot Expected value of CCCM and Variance
%
    clear all;
% 
    plot2D=0;
    col = ['r','b','g','y','m','r','b','g','y','m'];
    observable_name = 'M(dE_1,dE_2) CCCM: ';
    picsformat = 'png';
    workdir = pwd;
    ixdir = strfind(workdir,'/');
    dirname = workdir(ixdir(end-1)+1:ixdir(end)-1);
    emfig = ['EMdiag_',dirname];
    varfig= ['VarMdiag_',dirname];
    stdfig= ['StdMdiag_',dirname];
    
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
    M = zeros(p,p);
    V = zeros(p,p);
   
    ilev = 1;
    while ilev <= maxlevel
        fprintf(1,'    Level: %d \n',ilev);
        
        if ilev < maxlevel
            id = idxlevelE(ilev);
            file_name = file{id}.name;
            fM = fopen(file_name,'r');   % open the E(M)
            Ml = fread(fM,[p p],'double');                  % p x p matrix load
            ML{ilev}=Ml;
            fclose(fM);
            fprintf(1,'      File: %s \n',file_name);
            
            id = idxlevelDE(ilev);
            file_name = file{id}.name;
            fDM = fopen(file_name,'r');   % open the DE(M)
            DMl = fread(fDM,[p p],'double');                  % p x p matrix load
            fclose(fDM);
            fprintf(1,'      File: %s \n',file_name);
            
            id = idxlevelV(ilev);
            file_name = file{id}.name;
            fV  = fopen(file_name,'r');   % open the Var(M)
            M2l= fread(fV,[p,p],'double');
            VML{ilev} = M2l - Ml.^2;
            fclose(fV);
            fprintf(1,'      File: %s \n',file_name);
            
            id = idxlevelDV(ilev);
            file_name = file{id}.name;
            fDV  = fopen(file_name,'r');   % open the DVar(M)
            DM2l= fread(fDV,[p,p],'double');
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
            fM = fopen(file_name,'r');   % open the E(M)
            Ml = fread(fM,[p p],'double');                  % p x p matrix load
            ML{ilev}=Ml;
            fclose(fM);
            fprintf(1,'      File: %s \n',file_name);
            
            id = idxlevelV(ilev);
            file_name = file{id}.name;
            fV  = fopen(file_name,'r');   % open the Var(M)
            M2l= fread(fV,[p,p],'double');
            VML{ilev} = (M2l - Ml.^2)/(M_v(ilev)-1);
            fclose(fV);
            fprintf(1,'      File: %s \n',file_name);
            
            id = idxlevelENERGY(ilev);
            file_name = file{id}.name;
            fE = fopen(file_name,'r');   % open the file with energy mesh
            E = fread(fE,[p 1],'double');                   % p energy samps load
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
    % Plot the final value of E[M]
    strtitle = [observable_name, 'MLMC estimator: E_{MLMC}[M]'];
    ztitle = 'E[M(E_1,E_2)]';
    
    if plot2D == 1
        figure(1);clf;
        plot_CCCM(E,M,strtitle,ztitle,1);
    % 
    % Plot individual levels
        figure(2);clf;
        switch maxlevel
            case {2,3}
                for ilev = 1:maxlevel
                    strtitle = [observable_name, ' E_l[M] Level: l=',int2str(ilev)];
                    subplot(2,2,ilev);
                        plot_CCCM(E,ML{ilev},strtitle,ztitle,1);
                end 
            
                strtitle = [observable_name, ' E_{MLMC}[M]'];
                subplot(2,2,maxlevel+1);  
                plot_CCCM(E,M,strtitle,ztitle,1);
            otherwise
                spl = fix(sqrt(maxlevel+1));
                for ilev = 1:maxlevel
                    strtitle = [observable_name, ' E_l[M] Level: l=',int2str(ilev)];
                    subplot(spl+1,spl,ilev);
                     plot_CCCM(E,ML{ilev},strtitle,ztitle,1)
                end 
                strtitle = [observable_name, ' E_{MLMC}[M]'];
                subplot(spl+1,spl,maxlevel+1); 
                    plot_CCCM(E,M,strtitle,ztitle,1);
        end
    end
    
    ccmdiag = diag(M);
    for i=1:maxlevel
        Ml = ML{i};
        CCMdiag{i} = diag(Ml);
    end
    
    figure(3);clf;
    leg(1)=plot(E,ccmdiag/N_v(end)^2,'k','LineWidth',2); hold on;
    strleg{1} = 'E_{MLMC}[M] estimator';
    for ilev=1:maxlevel
        leg(ilev+1)=plot(E,CCMdiag{ilev}/N_v(ilev)^2,col(ilev)); hold on;
        strleg{ilev+1} = sprintf('E_l[M] l = %d',ilev);
    end
    xlabel('E'); ylabel('E[M(E,E)]');
    strtitle = [observable_name,' E[M(E,E)]']; 
    title(strtitle);
    legend(leg,strleg);
    annotation('textbox',[0.15, 0.7,0.3,0.18],...
        'String',{anns1,anns2,anns3,anns4},...
        'FontSize',9,'FitBoxToText','on');
    annotation('textbox',[0.15, 0.15,0.1,0.05],...
        'String',[dirname],...
        'FontSize',9,'LineStyle','none','FitBoxToText','on');
    saveas(gcf,emfig,picsformat);
   
    figure(4);clf;
       plot(ccmdiag/N_v(end)^2,E,'k','LineWidth',2); hold on;
       ylabel('E'); xlabel('E[M(E,E)]');
       strtitle = [observable_name,' E[M(E,E)] Poly order=',int2str(p)]; 
       title(strtitle);
       
 %----------------------------------------------------------------------

 %----------------------------------------------------------------------
 %
 % Plotting Variance estimates
 %
 % Plot the final value of E[M]
    strtitle = [observable_name, 'MLMC estimator: Var_{MLMC}[M]'];
    ztitle = 'Var[M(E_1,E_2)]';
    
    if plot2D == 1
        figure(5);clf;
        plot_CCCM(E,V,strtitle,ztitle,1);
    % 
    % Plot individual levels
        figure(6);clf;
        switch maxlevel
            case {2,3}
                for ilev = 1:maxlevel
                    strtitle = [observable_name, ' Var_l[M] Level: l=',int2str(ilev)];
                    subplot(2,2,ilev);
                        plot_CCCM(E,VML{ilev},strtitle,ztitle,1);
                end 
            
                strtitle = [observable_name, ' Var_{MLMC}[M]'];
                subplot(2,2,maxlevel+1);  
                plot_CCCM(E,V,strtitle,ztitle,1);
            otherwise
                spl = fix(sqrt(maxlevel+1));
                for ilev = 1:maxlevel
                    strtitle = [observable_name, ' Var_l[M] Level: l=',int2str(ilev)];
                    subplot(spl+1,spl,ilev);
                     plot_CCCM(E,VML{ilev},strtitle,ztitle,1)
                end 
                strtitle = [observable_name, ' Var_{MLMC}[M]'];
                subplot(spl+1,spl,maxlevel+1); 
                    plot_CCCM(E,V,strtitle,ztitle,1);
        end
    end
    
    varccmdiag = diag(V);
    for i=1:maxlevel
        VMl = VML{i};
        VarCCMdiag{i} = diag(VMl);
    end
%%    
    figure(7);clf;
    leg(1)=plot(E,varccmdiag,'k','LineWidth',2); hold on;
    strleg{1} = 'Var_{MLMC}[M] estimator';
    for ilev=1:maxlevel
        leg(ilev+1)=plot(E,VarCCMdiag{ilev},col(ilev)); hold on;
        strleg{ilev+1} = sprintf('Var_l[M] l = %d',ilev);
    end
    xlabel('E'); ylabel('Var[M(E,E)]');
    strtitle = [observable_name,' Var[M(E,E)]']; 
    title(strtitle);
    legend(leg,strleg);
    annotation('textbox',[0.15, 0.7,0.3,0.18],...
        'String',{anns1,anns2,anns3,anns4},...
        'FontSize',9,'FitBoxToText','on');
    annotation('textbox',[0.15, 0.15,0.1,0.05],...
        'String',[dirname],...
        'FontSize',9,'LineStyle','none','FitBoxToText','on');
     saveas(gcf,varfig,picsformat);
     
    figure(8);clf;
    leg(1)=semilogy(E,sqrt(varccmdiag),'k','LineWidth',2); hold on;
    strleg{1} = 'std_{MLMC}[M] estimator';
    for ilev=1:maxlevel
        leg(ilev+1)=semilogy(E,sqrt(VarCCMdiag{ilev}),col(ilev)); hold on;
        strleg{ilev+1} = sprintf('std_l[M] l = %d',ilev);
    end
    xlabel('E'); ylabel('std[M(E,E)]');
    strtitle = [observable_name,' std[M(E,E)]']; 
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
    