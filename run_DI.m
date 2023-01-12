clear all
%close all

I = 1;
eps_o = 0.9;
mu = 10;
theta_o = 0.5;
lambda = 10;

gamma = linspace(0, 5, 100);
gamma = gamma(2:end)';
C = linspace(0, I, 100);
C = C(2:end)';

Ui = zeros(length(gamma), length(C));
W = zeros(length(gamma), length(C));
Y = zeros(size(W));
Pay = zeros(size(W));

for i = 1:length(gamma)
    i
    for j = 1:length(C)
        [Ui(i, j), ~, ~, ~, ~, ~, W(i, j), Y(i, j), ~, ~, Pay(i, j), ~, ~, ~, ~, ~, ~, ~] = DI_exp(gamma(i), C(j), I, theta_o, eps_o,lambda, mu);
    end
end

fontsize = 20;

figure;
imagesc(C, gamma, W);
colormap(magma());
colorbar;
%caxis([0,1])
set(gca,'YDir','normal');
xlabel('$C$', 'interpreter', 'latex');
ylabel('$\gamma$', 'interpreter', 'latex');
title('Deterrence Effort ($W^*$)', 'interpreter', 'latex');
set(gca, 'FontSize', fontsize);
set(gcf,'Position', [100 100 600 500]);

figure;
imagesc(C, gamma, Y);
colormap(magma());
colorbar;
%caxis([0,1])
set(gca,'YDir','normal');
xlabel('$C$', 'interpreter', 'latex');
ylabel('$\gamma$', 'interpreter', 'latex');
title('Backup Effort ($Y^*$)', 'interpreter', 'latex');
set(gca, 'FontSize', fontsize);
set(gcf,'Position', [100 100 600 500]);

figure;
imagesc(C, gamma, Pay);
colormap(magma());
colorbar;
set(gca,'YDir','normal');
xlabel('$C$', 'interpreter', 'latex');
ylabel('$\gamma$', 'interpreter', 'latex');
title('Pay', 'interpreter', 'latex');
set(gca, 'FontSize', fontsize);
set(gcf,'Position', [100 100 600 500]);

figure;
imagesc(C, gamma, Ui);
colormap(magma());
colorbar;
set(gca,'YDir','normal');
xlabel('$C$', 'interpreter', 'latex');
ylabel('$\gamma$', 'interpreter', 'latex');
title('Ui', 'interpreter', 'latex');
set(gca, 'FontSize', fontsize);
set(gcf,'Position', [100 100 600 500]);