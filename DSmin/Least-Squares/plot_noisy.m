%Script for plotting integer least squares results.
%To be used for results with varying noise.


results_dir = '';

filefolderlist = dir(results_dir);
names = {filefolderlist.name};
name_mask = cellfun(@(x)contains(x, ".mat"), names);
file_names = natsortfiles(names(name_mask));

% fields = fieldnames(results);
snrs = [];

colors = distinguishable_colors(8); 
lines = ["-x", "--x", "-*", "--*", "-o", "-s", "-p", "-h", "-d", "-^", "-+", "-v"];

for i = 1:length(file_names)
    load(strcat(results_dir, "/", file_names{i}))
    snrs(i) = snr;
    
    ber_admm(i) = mean(1 - correct_admm);
    ber_rar(i) = mean(1 - correct_rar);
    ber_obq3(i) = mean(1 - correct_obq3);
    ber_ds(i) = mean(1 - correct_ds_rar);

    std_ber_admm(i) = std(1 - correct_admm);
    std_ber_rar(i) = std(1 - correct_rar);
    std_ber_obq3(i) = std(1 - correct_obq3);
    std_ber_ds(i) = std(1 - correct_ds_rar);

    rate_admm(i) = sum(correct_admm == 1) / n_trials;
    rate_rar(i) = sum(correct_rar == 1) / n_trials;
    rate_obq3(i) = sum(correct_obq3 == 1) / n_trials;
    rate_ds(i) = sum(correct_ds_rar == 1) / n_trials;

    if use_gurobi
        ber_gurobi(i) = mean(1 - correct_gurobi);
        std_ber_gurobi(i) = std(1 - correct_gurobi);
        rate_gurobi(i) = sum(correct_gurobi == 1) / n_trials;
    end


    optgap_admm(i) = mean((F_admm - Fopt) ./ Fopt);
    optgap_rar(i) = mean((F_rar - Fopt) ./ Fopt);
    optgap_obq3(i) = mean((F_obq3 - Fopt) ./ Fopt);
    optgap_ds(i) = mean((F_ds_rar - Fopt) ./ Fopt);

    std_optgap_admm(i) = std((F_admm - Fopt) ./ Fopt);
    std_optgap_rar(i) = std((F_rar - Fopt) ./ Fopt);
    std_optgap_obq3(i) = std((F_obq3 - Fopt) ./ Fopt);
    std_optgap_ds(i) = std((F_ds_rar - Fopt) ./ Fopt);

end

%%
clf;
ber_plot = figure;
set(0, 'CurrentFigure', ber_plot)
set(gca,'fontsize',20,'TickLabelInterpreter','latex')
hold on
std_ber_admm(2:2:length(std_ber_admm)) = NaN;
std_ber_rar(2:2:length(std_ber_rar)) = NaN;
std_ber_obq3(2:2:length(std_ber_obq3)) = NaN;
std_ber_ds(2:2:length(std_ber_ds)) = NaN;
errorbar(snrs, ber_admm, std_ber_admm, lines(1), 'DisplayName', "admm", 'color', colors(1, :),'linewidth', 2, 'markersize',10);
errorbar(snrs, ber_rar, std_ber_rar, lines(3), 'DisplayName', "rar", 'color', colors(2, :),'linewidth', 2, 'markersize',10);
errorbar(snrs, ber_obq3, std_ber_obq3, lines(5), 'DisplayName', "obq3", 'color', colors(3, :),'linewidth', 2, 'markersize',10);
errorbar(snrs, ber_ds, std_ber_ds, lines(6), 'DisplayName', "dsmin-rar", 'color', colors(4, :),'linewidth', 2, 'markersize',10);
if use_gurobi
    std_ber_gurobi(2:2:length(std_ber_gurobi)) = NaN;
    errorbar(snrs, ber_gurobi, [], std_ber_gurobi, lines(7), 'DisplayName', "Global", 'color', colors(5, :),'linewidth', 2, 'markersize',10);
end
% legend('Location','northeastoutside')
% legend(BackgroundAlpha=.7)
xlabel('$SNR_{dB}$', 'Interpreter','latex') 
ylabel('Bit Error Rate', 'Interpreter','latex') 
% ylim([0 max([ber_admm, ber_rar, ber_obq3, ber_ds])])
saveas(ber_plot,strcat(results_dir,'ber.fig'))
exportgraphics(gca,strcat(results_dir,'ber.pdf'))
hold off


clf;
rate_plot = figure;
set(0, 'CurrentFigure', rate_plot, 'defaulttextinterpreter','latex')
set(gca,'fontsize',20)
set(0,'DefaultTextFontname', 'Times')
set(0,'DefaultAxesFontName', 'Times')
% set(rate_plot, 'Units', 'pixels', 'Position', [10, 100, 1000, 400])
hold on
plot(snrs, rate_admm, lines(1), 'DisplayName', "ADMM", 'color', colors(1, :),'linewidth', 2, 'markersize',10);
plot(snrs, rate_rar, lines(3), 'DisplayName', "RAR", 'color', colors(2, :),'linewidth', 2, 'markersize',10);
plot(snrs, rate_obq3, lines(5), 'DisplayName', "OBQ", 'color', colors(3, :),'linewidth', 2, 'markersize',10);
plot(snrs, rate_ds, lines(6), 'DisplayName', "DCA", 'color', colors(4, :),'linewidth', 2, 'markersize',10);
if use_gurobi
    plot(snrs, rate_gurobi, lines(7), 'DisplayName', "Optimal", 'color', colors(5, :),'linewidth', 2, 'markersize',10);
end
legend('Location','southeast', BackgroundAlpha=.7)
xlabel('$SNR_{dB}$') 
ylabel('Recovery Probability') 
saveas(rate_plot,strcat(results_dir,'rate.fig'))
exportgraphics(gca,strcat(results_dir,'rate.pdf'))
hold off

clf;
f_plot = figure;
set(0, 'CurrentFigure', f_plot, 'defaulttextinterpreter','latex')
set(gca,'fontsize',20)
set(0,'DefaultTextFontname', 'Times')
set(0,'DefaultAxesFontName', 'Times')
hold on
std_optgap_admm(2:2:length(std_optgap_admm)) = NaN;
std_optgap_rar(2:2:length(std_optgap_rar)) = NaN;
std_optgap_obq3(2:2:length(std_optgap_obq3)) = NaN;
std_optgap_ds(2:2:length(std_optgap_ds)) = NaN;
errorbar(snrs, optgap_admm, std_optgap_admm, lines(1), 'DisplayName', "admm", 'color', colors(1, :),'linewidth', 2, 'markersize',10);
errorbar(snrs, optgap_rar, std_optgap_rar, lines(3), 'DisplayName', "rar", 'color', colors(2, :),'linewidth', 2, 'markersize',10);
errorbar(snrs, optgap_obq3, std_optgap_obq3, lines(5), 'DisplayName', "obq3", 'color', colors(3, :),'linewidth', 2, 'markersize',10);
errorbar(snrs, optgap_ds, std_optgap_ds, lines(6), 'DisplayName', "dsmin-rar", 'color', colors(4, :),'linewidth', 2, 'markersize',10);
% legend('Location','northeastoutside')
% legend(BackgroundAlpha=.7)
xlabel('$SNR_{dB}$') 
ylabel('$(F(\hat{x}) - F(x^{\natural}))/F(x^{\natural})$') 
saveas(f_plot,strcat(results_dir,'f.fig'))
exportgraphics(gca,strcat(results_dir,'f.pdf'))
hold off