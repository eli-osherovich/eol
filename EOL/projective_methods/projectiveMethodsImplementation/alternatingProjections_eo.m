classdef alternatingProjections_eo < hgsetget
    properties (SetAccess = protected)
        DampingFactor = 1;
        Projection1;
        Projection2;
    end
    
    % Optimization parameters (default value)
    properties 
        MaxIter = 100           % Maximal number of iterations (100)
        SaveDir = ''            % Directory to save iterations ('' = do not save)
        ComplexVarsFlag = false % Indicator whether the variable X is complex (false)
        Display = true          % Progress report (true)
    end
    
    
    methods
        function self = alternatingProjections_eo (proj1, proj2)
            self.Projection1 = proj1;
            self.Projection2 = proj2;
        end
        function x = run(self, x)
            
            % Save the original x's size.
            xSize = size(x);
            x = x(:);
            maxIter = self.MaxIter;
            saveDir = self.SaveDir;
            display = self.Display;
            proj1 = self.Projection1;
            proj2 = self.Projection2;
            
            % Create directory to save results (if required).
            if ~isempty(saveDir) && ~exist(saveDir, 'dir')
                mkdir(saveDir);
            end
            
            
            % Print progress report header (if required).
            if display
                fprintf('%10s %15s %15s %15s %15s\n', ...
                    'Iteration', 'Step Length','Total Err', 'Err 1', 'Err 2');
            end
            
            xOld = x;
            for i = 0:maxIter
                [x1, err1] = proj1.doProjection(x);
                [x2, err2] = proj2.doProjection(x1);
                
               
                
                % Save current x (if requested).
                if ~isempty(saveDir)
                    save(fullfile(saveDir, int2str(i)), 'x');
                end
                
                % Go to the new x.
                x = self.doStep(x, x1, x2);
                
                % Print progress (if requested).
                if display
                    fprintf('%10d %15.5e %15.5e %15.5e %15.5e\n', ...
                        i, norm(x - xOld), err1+err2, err1, err2);
                end
                
                xOld = x;
            end
            
            
            % Reshape x to its original size.
            x = reshape(x, xSize);
        end
    end
    
    methods (Access = protected)
        function x = doStep (self, x, ~, x2)
            % This class implements the simplest (damped) alternating
            % projections. Different subclasses will override this
            % function.
            
            df = self.DampingFactor;
            x = x + (x2-x) * df;
        end
    end
end

