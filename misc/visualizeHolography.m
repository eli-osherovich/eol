function visualizeHolography(xMag, refBeam, holImg, idx, phaseLo, phaseHi)

figure;
cla;
hold on;



% Origin.
O = complex(0, 0);

% Draw origin.
plot(O, 'ko', 'MarkerFaceColor', 'k');

% Draw the reference beam
plot([0; refBeam(idx)], 'k');
hold on;
axis square;

% Draw the circles corresponding to the reference beam (black), the
% measured intensity (red), and the unknown signal (green).
t = (0:0.001:1)*2*pi;
plot(abs(refBeam(idx))*cos(t), abs(refBeam(idx))*sin(t), 'k');
plot(sqrt(holImg(idx))*cos(t), sqrt(holImg(idx))*sin(t), 'r');
plot(real(refBeam(idx)) + ...
    xMag(idx)*cos(t), imag(refBeam(idx)) + xMag(idx)*sin(t), 'g');





yLo = sqrt(holImg(idx))*exp(1i*phaseLo(idx));
yHi = sqrt(holImg(idx))*exp(1i*phaseHi(idx));
y = (yLo+yHi)/2;
% plot segment
line_dir = (yHi - yLo) / abs(yHi - yLo);
hLo = plot([O;yLo], 'k');
hHi = plot([O;yHi], 'm');
legend([hLo, hHi], 'Low', 'High');

plot([O;y], ':');



