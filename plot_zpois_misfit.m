function plot_zpois_misfit(Z_saved, poisson_saved)
%h1 = histogram(Z_saved);
%h2 = histogram(poisson_saved);
h = histogram2(Z_saved,poisson_saved,100, 'displaystyle', 'Tile');
shading interp
%scatter3(Z_saved, poisson_saved, exp(log_posterior), 5, exp(log_posterior), 'o','filled');
%contour(Z_saved,poisson_saved,exp(log_posterior))
standard_figure
xlabel('Acoustic impedance')
ylabel('Poisson ratio')
colorbar