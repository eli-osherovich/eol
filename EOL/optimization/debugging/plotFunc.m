function  plotFunc(x0, dir, Ax0, Adir, FuncAxStruct, FuncXStruct, tInterval, nPoints)
% Plot function 

% Set default t Interval to [0, 1].
if nargin < 7
    tInterval = [0,1];
end

% Set default nPoints to 100.
if nargin < 8
    nPoints = 100;
end

% Space for Ax.
Ax = Ax0;

% Sample tInterval
allT = linspace(tInterval(1), tInterval(2), nPoints);

% Preallocate space for f.
f = zeros(numel(allT));

% Compute funciton values for all t values.
for i = 1:numel(allT)
    t = allT(i);
    
    x = x0 + t*dir;
    for k = 1:length(FuncAxStruct)
            Ax{k} = Ax0{k} + t * Adir{k};
    end
    f(i) = calc_EDx(x, Ax, FuncAxStruct, FuncXStruct, empty, [], false);
end

% Finally, plat f vs t.
figure(9999);
plot(allT, f);
grid on;
