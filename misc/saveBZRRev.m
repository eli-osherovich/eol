function saveBZRRev(saveFile, varargin)
    % Save BZR revision information
    
    
    % Copyright 2011 Eli Osherovich.
    
    
    if nargin < 2
        error('Not enough input arguments');
    end
    
    % Create empty cell for the revision numbers.
    revNumbers = cell(size(varargin));
    
    for i = 1:length(varargin)
        [status, output] = system(['bzr st ' varargin{i}]);
        
        % Check whether the command succeeded.
        if 0 ~= status 
            error('Failed to run BZR ST: %s', output);
        end
        
        % Non-empty OUTPUT means there are uncommited changes.
        if ~isempty(output)
            error('There are uncommited changes in %s', varargin{i});
        end
        
        % Save revision number of the current directory.
        [revStatus, revNum] = system(['bzr revno ' varargin{i}]);
        if 0 ~= revStatus
            error('Failed to run BZR REVNO: %s', revNum);
        end
        revNumbers{i} = revNum;
    end
    
    % Save all revisions (as a map container) to a file.
    revMap = containers.Map(varargin, revNumbers);
    save(saveFile, 'revMap');
    
