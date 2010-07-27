function val = getOpt_eo(options, opt_name, default_val)
% GETOPT_EO(OPTIONS, OPT_NAME, DEFAULT_VAL) - get option value.
% VAL = GETOPT_EO(OPTIONS, OPT_NAME, DEFAULT_VAL) - get option value named
% OPT_NAME from the struct OPTIONS. In case there is no such option use
% DEFAULT_VAL. The option name is case insensitive

% Examples:
%-------------------------------
% options.tolX = 0.0001;
% options.tolF = 0.1;
% options.Error = 'Error';
% options.error = 'error';
%
% GET AN EXISTING OPTION
% -----------------------------
% getOpt_eo(options, 'tolx', 1e-8)
% 
% ans =
% 
%    1.0000e-04
%
% GET A NON-EXISTING OPTION
% -----------------------------
% getOpt_eo(options, 'NE', 6)
% 
% ans =
% 
%      6
%
% TRYING TO GET AN OPTION THAT IS GIVEN TWICE
% --------------------------------------------
% getOpt_eo(options, 'error', 6)
% ??? Error using ==> getOpt_eo at 49
% Each option may appear only once : 'error' appears 2 times.


% Copyright 2010 Eli Osherovich.

% all given options
opt_names = fieldnames(options);

% find matches
matches = strcmpi(opt_name, opt_names);

switch sum(matches(:))
    
    case 0
        % not found: use default value
        val = default_val;
    case 1
        % one match found: use given value
        val = options.(opt_names{matches});
    otherwise
        % the option appears several times
        error('EOL:GETOPT:MultipleOptions',...
            'Each option may appear only once : ''%s'' appears %d times.',...
            opt_name, sum(matches(:)));
end
