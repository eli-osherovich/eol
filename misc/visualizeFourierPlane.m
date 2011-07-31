function visualizeFourierPlane(xTrueF, xRecF, idx, phaseLo, phaseHi)

figure;
cla;
hold on;

% True value's magnitude.
r = abs(xTrueF(idx));

% Origin.
O = complex(0, 0);

% Draw origin.
plot(O, 'ko', 'MarkerFaceColor', 'k');
% Draw circle.
t = (0:0.001:1)*2*pi;
plot(r*cos(t), r*sin(t), 'k');
axis square;

% Draw the original and reconstructed values
hXTrueF = plot(xTrueF(idx), 'or');
hXRecF = plot(xRecF(idx), 'sm');
legend([hXTrueF, hXRecF], 'True', 'Reconstructed');

yLo = r*exp(1i*phaseLo(idx));
yHi = r*exp(1i*phaseHi(idx));
y = (yLo+yHi)/2;
% plot segment
line_dir = (yHi - yLo) / abs(yHi - yLo);
hLo = plot([O;yLo], 'k');
hHi = plot([O;yHi], 'm');
legend([hLo, hHi], 'Low', 'High');

plot([O;y], ':');
plot([yLo; yLo + r * line_dir]);
plot([yHi; yHi - r * line_dir]);


