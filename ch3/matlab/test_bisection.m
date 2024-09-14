
clear a7  b7; clc; close all;
clf; hold off;

f = @(x) (sin(x) + sin(x/2))/4;

t = -5:0.01:9;
a0 = -5; b0 = 9;
ymin = -0.5 ; ymax = 0.5 ;

[a7, b7] = bisection_search(f, a0, b0, 5e-1);
rowp = size(a7,1);

figure(1);
%set(gcf,'Position',[0 0 900 400])
% Plot initial function 
axis([-5, 9, ymin, ymax]); grid;
hold on
    plot(t, f(t), 'LineWidth',2); 
    title("$f(x) = (\sin(x) + \sin(\frac{x}{2}))/4$", 'Interpreter','latex','FontSize',14);
    xlabel("$x$", 'Interpreter','latex','FontSize',14);
    ylabel("$f(x)$",'Interpreter','latex','FontSize',14);
    % The boundary is very narrow
    plot([a0 a0], [ymin, ymax], 'LineWidth',1,'Color','red'); 
    plot([b0 b0], [ymin, ymax], 'LineWidth',1, 'Color','red');
hold off

% Plot the bisection steps.
figure(2)
title("$f(x) = \frac{(\sin(x) + \sin(\frac{x}{2}))}{4}$", ...
    'Interpreter','latex','FontSize',14);

for i = 1: rowp
    subplot(1, rowp,i)
    axis([-5, 9, ymin, ymax]); grid;
    hold on
    plot(t, f(t), 'LineWidth',2); 
    xlabel("$x$", 'Interpreter','latex','FontSize',14);
    ylabel("$f(x)$",'Interpreter','latex','FontSize',14);
    % The boundary is very narrow
    plot([a7(i) a7(i)], [ymin, ymax], 'LineWidth',1,'Color','red'); 
    plot([b7(i) b7(i)], [ymin, ymax], 'LineWidth',1, 'Color','red');

    pgon1 = polyshape([ a7(i) a7(i) b7(i) b7(i) ], ...
        [ ymin, ymax, ymax, ymin ], 'Simplify',false);  
    plot(intersect(pgon1),'EdgeColor','none','FaceColor','b')
    

    hold off
end
drawnow;