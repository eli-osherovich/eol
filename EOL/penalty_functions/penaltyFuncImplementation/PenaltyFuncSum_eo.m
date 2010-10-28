classdef PenaltyFuncSum_eo < PenaltyFunc_eo
    
    properties (Access=private)
        PenaltyFuncListCell
    end
    
    methods 
        function self = PenaltyFuncSum_eo (varargin)
            % Make sure that all arguments are penalty functions
            if any(cellfun (@(c) ~isa(c, 'PenaltyFunc_eo'), varargin))
                error('EOL:PenaltyFuncSum:wrongArgType', ...
                    'Inputs must be of type PenaltyFunc_eo');
            end
            
            % Create a PenaltyFuncSum object
            self.PenaltyFuncListCell = varargin(:);
        end
   
        function [val, grad, hessMultVecorFunc] = doCalculations(self, x)
            pfList = self.PenaltyFuncListCell;
            
            nFunc = numel(pfList);
            
            val = 0;
            grad = 0;
            allHessMultFunc = cell(1, nFunc);
            allMultFactors = zeros(1, nFunc);
            
            % Since MATLAB does not call an overloaded SUBSREF inside
            % object's methods (and one of the entries in
            % PenaltyFuncListCell and be a PenaltySum_eo) we shall use an
            % explicit call. 
            
            % Prepare data for SUBSREF call
            subs.type = '()';
            subs.subs = {x};
            
            valCur = 0;
            gradCur = 0;
            for i = 1:nFunc
                pfCurr = pfList{i};
                                          
                switch nargout
                    case {0, 1}
                        valCur = subsref(pfCurr, subs);
                    case 2
                        [valCur, gradCur] = subsref(pfCurr, subs);
                    case 3
                        % Returning three outputs is a bit problematic.
                        % By the design, penalty fuctions should return (as
                        % the third output) a function that implements
                        % Hessian-vector product. However, this function
                        % knows nothing about the MULTFACTOR. Hence, we
                        % must take care of it.
                        
                        [valCur, gradCur, allHessMultFunc{i}] = doCalculations(pfCurr, x);
                        valCur = valCur * pfCurr.multFactor;
                        gradCur = gradCur * pfCurr.multFactor;
                        allMultFactors(i) = pfCurr.multFactor;
                end
                
                val = val + valCur;
                grad = grad + gradCur;
            end
            
            hessMultVecorFunc  = @hessMult;
            
            function hessV = hessMult(v)
                hessV = 0;
                for k = 1:nFunc
                    hessV = hessV + allMultFactors(k)*allHessMultFunc{k}(v);
                end
           end
        end
    end
end
