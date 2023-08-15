function turbineMap(Tt4, Pt4, N, geom)

Mdat = load("Mdat.mat");

%% Turbine Inlet Conditions

mdotc = 10:0.1:60;
Nc = 1700:500:3200;
Nc(end + 1) = N/sqrt(Tt4/288.15);

% Determine efficienct and pressure ratio for ranges of corrected mass flow
% rate and N.
for j = 1:length(Nc)
    for i=1:length(mdotc)
        [eta_t(i,j), Pr(i,j)] = turb(mdotc(i), Nc(j), Tt4, Pt4, geom, Mdat);
        Pr(i,j) = 1/Pr(i,j);
    end
end

%% Plots
figure (10)
hold on
for j = 1:length(Nc) - 1
    plot(mdotc, Pr(:,j), 'LineWidth', 1);
    leg{j} = sprintf('$N/\\sqrt{\\theta}$=%d', Nc(j));
end
plot(mdotc, Pr(:,end), 'LineWidth', 2, 'LineStyle','--', 'Color','black');
leg{j + 1} = sprintf('$N/\\sqrt{\\theta}$=%d', floor(Nc(end)));
legend(leg, 'Interpreter','latex', 'Location','best');
title('Turbine Performance Map')
xlabel('$\dot{m}\frac{\sqrt{\theta_4}}{\delta_4}$', 'Interpreter','latex')
ylabel('$\frac{P_{t4}}{P_{t5}}$', 'Interpreter','latex')
fontsize(gca,14,"points")
set(gcf, 'Position',  [500, 200, 800, 600])
hold off

figure (11)
hold on
for j = 1:length(Nc) - 1
    plot(mdotc, eta_t(:,j), 'LineWidth', 1);
    leg{j} = sprintf('$N/\\sqrt{\\theta}$=%d', Nc(j));
end
plot(mdotc, eta_t(:,end), 'LineWidth', 2, 'LineStyle','--', 'Color','black');
leg{j + 1} = sprintf('$N/\\sqrt{\\theta}$=%d', floor(Nc(end)));
legend(leg, 'Interpreter','latex', 'Location','best');
title('Turbine Performance Map')
xlabel('$\dot{m}\frac{\sqrt{\theta_4}}{\delta_4}$', 'Interpreter','latex')
ylabel('$\eta_t$', 'Interpreter','latex')
fontsize(gca,14,"points")
set(gcf, 'Position',  [500, 200, 800, 600])
hold off