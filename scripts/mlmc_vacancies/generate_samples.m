function generate_samples(pVacs,biLev,N1,N2,nrat,orbindcs,nrsamples,dirstr)
% function GENERATE_SAMPLES samples vacancies for one level at a time
% Samples are written to text files in the format required by LEO2D with
% explicit diagonalization. That is, for each sample a "JOB=nr:" line 
% followed by a line with a comma separated list of integers denoting the 
% orbitals that are removed by the vacancies, and then an empty line before
% the next sample.
%
% INPUT:
%   pVacs: vector of reals on [0,1] given vacancy probability per atom site
%   biLev: true if this is a bi-level, i.e. create control variate samples
%   N1,N2: Positive integers for super cell size in either direction
%   nrat:  Positive integer for nr of atoms per unit cell 
%   orbindcs: Vector with values in {1,...,nrat} indicating which atom an
%             orbital belongs to
%   nrsamples: Positive integer for number of random samples
%   dirstr: String with name of directory to save output
%
% No output - results written to file

  check_input(pVacs,biLev,N1,N2,nrat,orbindcs)
  folder_name = make_folder(dirstr);
  
  filename = 'vacancies.dat';
  vacancy_string = [];

  if biLev
    filenames_CV = repmat('vacancies_CVX.dat',4,1);
    for k=1:4,
      strrep(filenames_CV(k,:),'X',int2str(k));
    end
    vacancy_CV_strings = {[],[],[],[]};
  end
  for k=1:nrsamples
    [atomvac,orbvac] = generate_Vacancies(pVacs,N1,N2,nrat,orbindcs);
    vacancy_string = [vacancy_string,Vac_to_str(orbvac,k),'\n'];
    if biLev % Create vacancy lists for the control variates
      orbvacs_CV = generate_CV_vacancies(N1,N2,nrat,orbindcs,atomvac);
      for l=1:4
        vacancy_CV_strings{l} = [vacancy_CV_strings{l},...
          Vac_to_str(orbvacs_CV{l},k),'\n'];
      end
    end
  end
  fh = fopen([folder_name,filesep,filename],'w+t');
  fprintf(fh,vacancy_string);
  fclose(fh);
  if biLev
    filename_CV = strsplit(filename,'.');
    for l=1:4
      fh = fopen([folder_name,filesep,filename_CV{1},'_CV',int2str(l),'.',filename_CV{2}],'w+t');
      fprintf(fh,vacancy_CV_strings{l});
      fclose(fh);
    end    
  end
end

function check_input(pVacs,biLev,N1,N2,nrat,orbindcs)
  if biLev && (rem(N1,2) || rem(N2,2))
    error('N1 and N2 must be divisible by 2 in bi-level case')
  end
  if biLev && N1~=N2
    error('Multilevel generation only implemented for N1=N2')
  end
  if min(pVacs)<0 || max(pVacs)>1
    error('Bad Probabilities')
  end
  if length(unique(orbindcs))~=nrat || max(orbindcs)>nrat || min(orbindcs)<1
    error('Range of orbital indices must equal 1,...,nrat')
  end
end

function [atomvac,orbvac] = generate_Vacancies(pVacs,N1,N2,nrat,orbindcs)
  atomvac = rand(nrat,N1*N2)<repmat(pVacs,1,N1*N2); % true = remove atom
  nrorb   = length(orbindcs); % orbitals per unit cell
  orbvac = [];
  for k=1:nrat
    myorbs = find(orbindcs==k);
    myvac  = find(atomvac(k,:));
    orbvac = sort([orbvac,repmat(myorbs,1,length(myvac)) + ...
      kron((myvac-1)*nrorb,ones(1,length(myorbs)))]); % a bit ugly...
  end
end

function orbvacs_CV = generate_CV_vacancies(N1,N2,nrat,orbindcs,atomvac)
  % divide into four corners
  icell = 1:N1*N2;
  iCV1 = mod(icell-1,N1)<N1/2 & (icell-1)/N1<N2/2;
  iCV2 = ~iCV1 & (icell-1)/N1<N2/2;
  iCV3 = ~iCV1 & mod(icell-1,N1)<N1/2;
  iCV4 = ~(iCV1|iCV2|iCV3);
  atomvac = [atomvac(:,iCV1);atomvac(:,iCV2);atomvac(:,iCV3);atomvac(:,iCV4)];

  nrorb   = length(orbindcs); % orbitals per unit cell
  orbvacs_CV = cell(1,4);
  for l=1:4
    orbvac = [];
    for k=1:nrat
      myorbs = find(orbindcs==k);
      myvac  = find(atomvac((l-1)*nrat+k,:));
      orbvac = sort([orbvac,repmat(myorbs,1,length(myvac)) + ...
        kron((myvac-1)*nrorb,ones(1,length(myorbs)))]);
    end
    orbvacs_CV{l} = orbvac;
  end
  
end

function outstr = Vac_to_str(orbvac,nr)
  % ugly
  outstr = ['JOBID = ',int2str(nr),'\nCLUSTERID = -1\nVAC = '];
  if isempty(orbvac) == 1
    outstr = [outstr,'-1'];  
  else
    for k=1:length(orbvac), 
      outstr = [outstr,int2str(orbvac(k)),','];
    end
  end
  if outstr(end)==',', 
    outstr(end) = [];
  end
  outstr = strcat(outstr,'\nTAR = -1\n');
end

function target_dir = make_folder(dirstr)
  outdir_base = 'Vacancies';
  
  if exist(outdir_base,'file') && ~exist('Vacancies','dir') % dir is a file
    loc_dir = dir;
    for k=1:length(loc_dir)
      if strcmp(loc_dir(k).name,outdir_base)
        error(['''',outdir_base,...
          ''' exists in current directory and is not a directory'])
      end
    end
  elseif ~exist(outdir_base,'dir')
    eval(['!mkdir ',outdir_base]);
  end
  
  if isempty(dirstr)
    ctime = clock;
    target_dir = [outdir_base,filesep,'Vac_',int2str(ctime(1)),'_',...
      int2str(ctime(2)),'_',int2str(ctime(3)),'__',int2str(ctime(4)),...
      '_',int2str(ctime(5))];
  else
    target_dir = [outdir_base,filesep,'Vac_',dirstr];
  end

  if exist(target_dir,'dir')
    overwrite = input('Target directory already exists. Overwrite? (y/n)','s');
    if ~strcmp(overwrite,'y')
      error('Exited on user request')
    end
  else
    eval(['!mkdir ',target_dir]);
  end
  
  % Since "make_folder" is called before the random vacancies are sampled,
  % saving the state of the random number generator allows for systematic
  % debugging among other things.
  randnstate = rng;
  save([target_dir,filesep,'RNG_initial_state'],'randnstate')
end