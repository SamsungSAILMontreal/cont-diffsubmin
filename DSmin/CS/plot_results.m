results_dir = '/Users/georgeorfanides/Desktop/thesis/DSLattice/CS/omp/signal_length=256,sparsity=26,n_trials=100,jobid=123/';

filefolderlist = dir(results_dir);
filteredlist = filefolderlist(~ismember({filefolderlist.name}, {'.', '..', 'l2error.fig', 'suppdiff.fig', 'success.fig', '.DS_Store'}));

file_names = {filteredlist.name};

L = regexp(results_dir,'(?<=signal_length=)\d*','match');

signal_length = str2num(L{1});

colors = distinguishable_colors(length(file_names)); 
% colors = [colors(1, :); colors(1, :);colors(2, :); colors(2, :); colors(3:end, :)];
lines = ["-x", "--x", "-*", "--*", "-o", "-s", "-p", "-h", "-d", "-^", "-+", "-v"];
names = ["OMP", "BOX-LASSO", "DCA", "DCA-LASSO"];

d = dictionary("omp", "OMP",...
                "lasso-roundlattice-warmstart", "LASSO",...
                "dsminLocalSearch-warmstart", "DCA",...
                "dsminLocalSearch-warmstart-LassoInit", "DCA-LASSO", ...
                "dsminLocalSearch-warmstart-randomperm", "DCA-Random");

color_dict = dictionary("omp", 	"#0000FF",...
                "lasso-roundlattice-warmstart", "#FF0000",...
                "dsminLocalSearch-warmstart", "#000000",...
                "dsminLocalSearch-warmstart-LassoInit", "#000000", ...
                "dsminLocalSearch-warmstart-randomperm", "#FF0000");

line_dict = dictionary("omp", "-x",...
                "lasso-roundlattice-warmstart", "-*",...
                "dsminLocalSearch-warmstart", "-o",...
                "dsminLocalSearch-warmstart-LassoInit", "--o", ...
                "dsminLocalSearch-warmstart-randomperm", "-s");

eps = 1e-5;

clf;
l2_plot = figure;
hold(gca,'on')

supp_plot = figure;
hold(gca,'on')

success_plot = figure;
hold(gca,'on')
for k = 1:length(file_names)
    load(strcat(results_dir, file_names{k}))
    method = results.method;
    field_names = string(fieldnames(results));
    field_mask = cellfun(@(x)contains(x, "n_obs"), fieldnames(results));
    result_fields = field_names(field_mask);
    n_obs = [];
    for i = 1:length(result_fields)
        n_obs(i) = str2num(erase(result_fields(i), 'n_obs'));
    end

    len = length(n_obs);
    
    l2error_results = [];
    support_results = [];
    success_results = [];
    std_l2error = [];
    std_support = [];
    for i = 1:length(result_fields)
        l2error = results.(result_fields(i)).l2error;
        support_diff = results.(result_fields(i)).suppdiff;
        N = length(l2error);
        l2error_results(i) = mean(l2error);
        support_results(i) = mean(support_diff);
        success_results(i) = mean(l2error <= 1e-3);
        std_l2error(i) = std(l2error);
        std_support(i) = std(support_diff);
    end
    set(0, 'CurrentFigure', l2_plot)
%     plot(n_obs / signal_length, l2error_results, lines(k), 'DisplayName', method, ...
%         'color', colors(k, :),'linewidth', 2, 'markersize',10);
    % std_l2error(2:2:length(std_l2error)) = NaN;
    upper_std_l2error = nan(1,len);
    upper_std_l2error(2:2:len) = std_l2error(2:2:len);
    lower_std_l2error = nan(1,len);
    lower_std_l2error(2:2:len) = min(l2error_results(2:2:len), std_l2error(2:2:len));
    errorbar(n_obs / signal_length, l2error_results, lower_std_l2error, upper_std_l2error, line_dict(method), 'DisplayName', d(method), ...
        'color', color_dict(method),'linewidth', 2, 'markersize',10);

    set(0, 'CurrentFigure', supp_plot)
%     plot(n_obs / signal_length, support_results, lines(k), 'DisplayName', method, ...
%         'color', colors(k, :),'linewidth', 2, 'markersize',10);
    % std_support(2:2:length(std_support)) = NaN;
    upper_std_support = nan(1,len);
    upper_std_support(2:2:len) = std_support(2:2:len);
    lower_std_support = nan(1,len);
    lower_std_support(2:2:len) = min(support_results(2:2:len), std_support(2:2:len));
    errorbar(n_obs / signal_length, support_results, lower_std_support, upper_std_support, line_dict(method), 'DisplayName', d(method), ...
            'color', color_dict(method),'linewidth', 2, 'markersize',10);

    set(0, 'CurrentFigure', success_plot)
%     plot(n_obs / signal_length, support_results, lines(k), 'DisplayName', method, ...
%         'color', colors(k, :),'linewidth', 2, 'markersize',10);
    plot(n_obs / signal_length, success_results, line_dict(method), 'DisplayName', d(method), ...
            'color', color_dict(method),'linewidth', 2, 'markersize',10);
end

figure(l2_plot)
set(gca,'fontsize',20,'TickLabelInterpreter','latex')
% set(l2_plot, 'Units', 'pixels', 'Position', [10, 100, 1000, 400])
% legend('Location','northeastoutside')
xlabel('m / n' , 'Interpreter','latex') 
ylabel('Estimation Error', 'Interpreter','latex') 
xlim([0.1 1])
saveas(l2_plot,strcat(results_dir,'l2error.fig'))
exportgraphics(gca,strcat(results_dir,'l2error.pdf'))

figure(supp_plot)
set(gca,'fontsize',20,'TickLabelInterpreter','latex')
% set(supp_plot, 'Units', 'pixels', 'Position', [10, 100, 1000, 400])
% legend('Location','northeastoutside')
xlabel('m / n', 'Interpreter','latex') 
ylabel('Support Error', 'Interpreter','latex') 
xlim([0.1 1])
saveas(supp_plot,strcat(results_dir,'suppdiff.fig'))
exportgraphics(gca,strcat(results_dir,'suppdiff.pdf'))

figure(success_plot)
set(gca,'fontsize',20,'TickLabelInterpreter','latex')
% set(success_plot, 'Units', 'pixels', 'Position', [10, 100, 1000, 400])
legend('Location','northwest', BackgroundAlpha=.7)
xlabel('m / n', 'Interpreter','latex') 
ylabel('Recovery Probability', 'Interpreter','latex') 
xlim([0.1 1])
saveas(success_plot,strcat(results_dir,'success.fig'))
exportgraphics(gca,strcat(results_dir,'success.pdf'))