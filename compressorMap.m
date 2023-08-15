function compressorMap(Tt2, Pt2, N, geom)

Mdat = load("Mdat.mat");

mdotc = 100:300;
Nc = 2000:500:6500;
Nc(end + 1) = N/sqrt(Tt2/288.15);
Nc = unique(Nc);
% Index to find operating N
indsp = find(Nc == N, 1);

% Determine efficienct and pressure ratio for ranges of corrected mass flow
% rate and N.
for j = 1:length(Nc)
    for i=1:length(mdotc)
%         try
        [eta_c(i,j), Pr(i,j), ~, ~, ~, ~] = compr(mdotc(i), Nc(j), Tt2, Pt2, geom, Mdat);
%         catch
%             disp('here')
%         end
    end
end

%% Contours
Ncex = Nc(1):Nc(end);
[A, B] = meshgrid(Ncex, mdotc);
len1 = size(A,1);
len2 = size(A,2);
parfor i=1:len1
    for j = 1:len2
        [coneta_c(i,j), conPr(i,j)] = compr(B(i,j), A(i,j), Tt2, Pt2, geom, Mdat);
    end
end

%% Plotting Specifics

% Values and indices of stall line, and definitions of annotations
[stall, ind] = max(Pr);
mdotcStall = mdotc(ind);

% Values and indices of choke line, and definitions of annotations
[row, col] = find(~isnan(Pr));
mdotcChoke = mdotc(row(end));
Prchoke = Pr(row(end), col(end));


%% Plots
figure (1)
hold on
for j = 1:length(Nc)
    if j == indsp
        plot(mdotc, Pr(:,indsp), 'LineWidth', 2, 'LineStyle','--', 'Color','blue');
        leg{j} = sprintf('$N/\\sqrt{\\theta}$=%d', floor(Nc(indsp)));
        continue;
    end
    plot(mdotc, Pr(:,j), 'LineWidth', 1);
    leg{j} = sprintf('$N/\\sqrt{\\theta}$=%d', Nc(j));
end
leg{end+1} = sprintf('Surge Line');
leg{end+1} = sprintf('Choke Line');
leg{end+1} = sprintf('Windmilling Line');
plot(mdotcStall, stall, 'LineWidth', 2, 'Marker','o', 'Color','black');
plot(mdotcChoke*ones(10, 1), linspace(1, Prchoke, 10), 'LineWidth', 2, 'Marker','x', 'Color','black');
plot(linspace(mdotc(row(1)), mdotc(row(end)), 10), ones(10,1), 'LineWidth', 2, 'Marker','>', 'Color','black');
[C,h] = contour(B,conPr, coneta_c, 'LineWidth', 0.25, 'LineColor','black');
clabel(C,h);
legend(leg, 'Interpreter','latex', 'Location','best', 'FontSize', 10);
title('Compressor Performance Map')
xlabel('$\dot{m}\frac{\sqrt{\theta_2}}{\delta_2}$', 'Interpreter','latex')
ylabel('$\frac{P_{t3}}{P_{t2}}$', 'Interpreter','latex')
ylim([0.5 20])
fontsize(gca,14,"points")
set(gcf, 'Position',  [500, 200, 800, 600])
hold off