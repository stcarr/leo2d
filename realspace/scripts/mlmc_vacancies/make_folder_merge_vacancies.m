function target_dir = make_folder_merge_vacancies
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
  
  ctime = clock;
  target_dir = [outdir_base,filesep,'MLMCRun_',int2str(ctime(1)),'_',...
    int2str(ctime(2)),'_',int2str(ctime(3)),'_',int2str(ctime(4)),...
    '_',int2str(ctime(5))];
  if exist(target_dir,'dir')
    overwrite = input('Target directory already exists. Overwrite? (y/n)','s');
    if ~strcmp(overwrite,'y')
      error('Exited on user request')
    end
  else
    eval(['!mkdir ',target_dir]);
  end
  
end