function saveBZRRev(saveFile, varargin)
    % Save BZR revision information.
    % SAVEBZRREV(SAVEFILE) - save BZR revision information of the calling
    % function and (all) its callers.
    % SAVEBZRREV(SAVEFILE, DIR1, DIR2, ...) in addition to the callers of
    % the current function include directories DIR1, DIR2, etc. to the
    % report.
    % SAVEBZRREV('',...) - if SAVEFILE is empty the function will not save
    % any information all test, however will run as usual.
    
    
    
    % Copyright 2011 Eli Osherovich.
    
    
    % Get all callers.
    callerStruct = dbstack(1, '-completenames');
    callerDirs = cellfun(@fileparts, ... 
        {callerStruct.file}, 'UniformOutput', false);
    
    % Create empty cell for the revision numbers.
    revNumbers = cell(numel(callerDirs) +  numel(varargin), 1);
    
    % All dirs.
    allDirs = [callerDirs, varargin];
    
    for i = 1:length(allDirs)
        [status, output] = system(['bzr st ' allDirs{i}]);
        
        % Check whether the command succeeded.
        if 0 ~= status 
            error('Failed to run BZR ST: %s', output);
        end
        
        % Non-empty OUTPUT means there are uncommited changes.
        if ~isempty(output)
            error('There are uncommited changes in %s', allDirs{i});
        end
        
        % Save revision number of the current directory.
        [revStatus, revNum] = system(['bzr revno ' allDirs{i}]);
        if 0 ~= revStatus
            error('Failed to run BZR REVNO: %s', revNum);
        end
        revNumbers{i} = revNum;
    end
    
    % Save all revisions (as a map container) to a file.
    if ~isempty(saveFile)
        revMap = containers.Map(allDirs, revNumbers);
        save(saveFile, 'revMap');
    end
    
