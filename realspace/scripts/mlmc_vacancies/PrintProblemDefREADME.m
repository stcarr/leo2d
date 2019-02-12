%
% Problem definition
%
% Print to a file defined as freadme
%
clear all;
%
freadme = fopen('../../README','a');
% 
    if exist('./problem_definition.mat','file') == 2
        %
        % Config file stored during generating (newer versions of problem
        % definitions)
        %
        load problem_definition.mat;
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
        
%
% Store info to a file
        dirname = pwd;
        idir = strfind(dirname,'/');
        dname = dirname(idir(end-1)+1:idir(end)-1);
        
        fprintf(freadme,'%s    Number of levels:      %d  of %d\n',dname,nlev,Nr_ext+1);
        fprintf(freadme,'   Polynomial order:      %d \n', p );
        fprintf(freadme,'   Defect probabilities:  Mo: %4.3f  X_A: %4.3f  X_B: %4.3f \n',pVacs);
        fprintf(freadme,'   Levels:\n');
        for i=1:nlev
            fprintf(freadme,'   l = %d size = %d x %d  # samples = %d \n',i,...
                N_v(i),N_v(i),nslev(i));
        end
        fprintf(freadme,'\n\n');
    else
        % 
        % Config file not stored (older versions)
        % Set the values manually
        fprintf(1,'WARNING: problem_definition.mat not found \n');
        fprintf(1,'         check that the set values are correct\n');
    end
    
    fclose(freadme);
    