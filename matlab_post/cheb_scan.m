  function [input_data, data, dos, samples] = cheb_scan()

    workspace;  % Make sure the workspace panel is showing.
    format longg;
    format compact;
    
    nShifts = 1;
    poly_order = 4000;
    energy_rescale = 20;

    % Define a starting folder.
    start_path = fullfile('C:\Users\Stephen\Documents\MobaXterm\home');
    % Ask user to confirm or change.
    topLevelFolder = uigetdir(start_path);
    
    %jobPool = parpool(4);

    if topLevelFolder == 0
        return;
    end
    % Get list of all subfolders.
    allSubFolders = genpath(topLevelFolder);
    % Parse into a cell array.
    remain = allSubFolders;
    listOfFolderNames = {};
    while true
        [singleSubFolder, remain] = strtok(remain, ';');
        if isempty(singleSubFolder)
            break;
        end
        listOfFolderNames = [listOfFolderNames singleSubFolder];
    end
    numberOfFolders = length(listOfFolderNames);
    outCount = 0;
    find_targets = 0;

    % Process all image files in those folders.
    for k = 1 : numberOfFolders
        % Get this folder and print it out.
        thisFolder = listOfFolderNames{k};
        fprintf('Processing folder %s\n', thisFolder);

        % Get PNG files.
        filePattern = sprintf('%s/*.cheb', thisFolder);
        baseFileNames = dir(filePattern);
        numberOfOutFiles = length(baseFileNames);
        % Now we have a list of all files in this folder.

        if numberOfOutFiles >= 1
            
            inputID = fopen(sprintf('%s/hstruct.in', thisFolder));
            inputStrings = textread(sprintf('%s/hstruct.in', thisFolder),'%s');
            for i = 1:length(inputStrings)
                if strcmp('NUM_SHEETS',inputStrings(i)) == 1
                    num_sheets = str2num(cell2mat(inputStrings(i+2)));
                    unitCells = zeros(3,3,num_sheets);
                    alphas = zeros(num_sheets,1);
                end
                    
                if strcmp('START_SHEET',inputStrings(i)) == 1
                    current_sheet = str2num(cell2mat(inputStrings(i+1)));
                end
                if strcmp('ALPHA',inputStrings(i)) == 1
                    alpha(current_sheet) = str2num(cell2mat(inputStrings(i+2)));
                end
                
                if strcmp('UNITCELL1',inputStrings(i)) == 1
                    unitCells(1,1,current_sheet) = str2num(cell2mat(inputStrings(i+2)));
                    unitCells(1,2,current_sheet) = str2num(cell2mat(inputStrings(i+3)));
                    unitCells(1,3,current_sheet) = str2num(cell2mat(inputStrings(i+4)));
                end
                
                if strcmp('UNITCELL2',inputStrings(i)) == 1
                    unitCells(2,1,current_sheet) = str2num(cell2mat(inputStrings(i+2)));
                    unitCells(2,2,current_sheet) = str2num(cell2mat(inputStrings(i+3)));
                    unitCells(2,3,current_sheet) = str2num(cell2mat(inputStrings(i+4)));
                end
                
                if strcmp('UNITCELL3',inputStrings(i)) == 1
                    unitCells(3,1,current_sheet) = str2num(cell2mat(inputStrings(i+2)));
                    unitCells(3,2,current_sheet) = str2num(cell2mat(inputStrings(i+3)));
                    unitCells(3,3,current_sheet) = str2num(cell2mat(inputStrings(i+4)));
                end
                
                if strcmp('NUM_ORBITALS',inputStrings(i)) == 1
                    num_orbitals(current_sheet) = str2num(cell2mat(inputStrings(i+2)));
                end
                
                if strcmp('POS',inputStrings(i)) == 1
                    for j = 1:num_orbitals
                        for k = 1:3
                            m = i + 1 + (j-1)*3 + k;
                            pos(j,k,current_sheet) = str2num(cell2mat(inputStrings(m)));
                        end
                    end
                end
                
                if strcmp('HEIGHT',inputStrings(i)) == 1
                    heights(current_sheet) = str2num(cell2mat(inputStrings(i+2)));
                end
                
                if strcmp('ANGLE',inputStrings(i)) == 1
                    angles(current_sheet) = str2num(cell2mat(inputStrings(i+2)));
                end
                
                if strcmp('POLY_ORDER',inputStrings(i)) == 1
                    poly_order = str2num(cell2mat(inputStrings(i+2)));
                end
                
                if strcmp('ENERGY_RESCALE',inputStrings(i)) == 1
                    energy_rescale = str2num(cell2mat(inputStrings(i+2)));
                end
                
                if strcmp('NSHIFTS',inputStrings(i)) == 1
                    nShifts = str2num(cell2mat(inputStrings(i+2)));
                end
            end
            
            fclose(inputID);
            
            inputData(1,1) = nShifts;
            inputData(2,1) = poly_order;
            inputData(3,1) = num_sheets;
            
            data_indexing = 4;
            for i = 1:num_sheets;
                inputData(data_indexing:data_indexing+2,1:3) = unitCells(:,:,i);
                inputData(data_indexing+3,1) = heights(i);
                inputData(data_indexing+4,1) = angles(i);
                inputData(data_indexing+5,1) = alpha(i);
                inputData(data_indexing+6,1) = num_orbitals(i);
                data_indexing = data_indexing + 7;
                for j = 1:num_orbitals(i)
                    inputData(data_indexing,1:3) = pos(j,1:3,i); 
                    data_indexing = data_indexing + 1;
                end
                
            end
            
            raw_avg_dos = 0;
            % Go through all files.
            for f = 1 : numberOfOutFiles
                
                fullFileName = fullfile(thisFolder, baseFileNames(f).name);
                fprintf('     Processing file %s\n', fullFileName);
                in_file = fopen(fullFileName);
                
                in_line = fgetl(in_file);
                
                while ischar(in_line)
                    
                    inputStrings = strsplit(in_line,{' ',','},'CollapseDelimiters',true);
                    
                    if (isempty(inputStrings) == 1)
                        continue
                        
                    end
                    %fprintf('%s \n',inputStrings{1}); 

                    if strcmp('JOBID',inputStrings{1}) == 1
                        
                        temp_data_index = data_indexing;
                        jobID = str2num(inputStrings{1+2});
                        fprintf('Working on JOBID = %d\n', jobID);
                        inputData(temp_data_index,1) = jobID;
                        temp_data_index = temp_data_index + 1;
                        outCount = outCount + 1;
                        find_targets = 0;
                        find_shifts = 0;

                    end
                    
                    if ( strcmp('SHEET:',inputStrings{1}) == 1)
                        find_shifts = 1;
                    end
                    
                    if (find_shifts == 1)
                        for s = 1:num_sheets;
                            if ( strcmp(num2str(s),inputStrings{1}) )
                                shifts(s,1) = str2num(inputStrings{3});
                                shifts(s,2) = str2num(inputStrings{4});
                                shifts(s,3) = str2num(inputStrings{5});
                                
                                if (s == num_sheets)
                                    find_shifts = 0;
                                    inputData(temp_data_index:temp_data_index+(num_sheets-1),1:3) = shifts(:,:);
                                    temp_data_index = temp_data_index + num_sheets;
                                end
                            end
                        end
                    end
                    
                    if strcmp('NUM_TAR',inputStrings{1}) == 1

                        num_tar = str2num(inputStrings{1+2});
                        inputData(temp_data_index,1) = num_tar;
                        temp_data_index = temp_data_index + 1;

                    end
                    if strcmp('TAR_LIST:',inputStrings{1}) == 1
                        
                        tar_list = zeros(num_tar,1);
                        for tar = 1:num_tar
                           tar_list(tar) =  str2num(inputStrings{1+tar});
                           inputData(temp_data_index,tar) = tar_list(tar);
                        end
                        
                        temp_data_index = temp_data_index + 1;
                    end
                    
                    if ( strcmp('T:',inputStrings{1}) == 1)
                        find_targets = 1;
                    end

                    if (find_targets == 1)
                        for x = 1:size(tar_list,1)

                            if ( strcmp(strcat(num2str(tar_list(x)),':'), inputStrings{1}) == 1 )

                                raw_data = zeros(poly_order,1);
                                for d_index = 1:(poly_order-1)
                                    raw_data(d_index) = str2num(inputStrings{1+d_index});
                                end

                                [t_d,t_s] = cheb_compute(raw_data,1,poly_order,energy_rescale);
                                temp_dos(x,:) = t_d;
                                temp_samples(x,:) = t_s;
                                
                                if (x == size(tar_list,1))
                                
                                    raw_dos(outCount) = cellCreate(temp_dos);
                                    raw_samples(outCount) = cellCreate(temp_samples);
                                    raw_input_data(outCount) = cellCreate(inputData);

                                end
                            end
                          
                            
                        end
                    end
                    
                    in_line = fgetl(in_file);
                end
            end
        end
    end

    data = raw_data;
    dos = raw_dos;
    input_data = raw_input_data;
    samples = raw_samples;
    
    %delete(jobPool);
end

function [cellOut] = cellCreate(mat_in)
    [r,c] = size(mat_in);
    cellOut = mat2cell(mat_in,r,c);
end