function plot_zbeta_posterior(Z_saved, beta_saved, log_posterior_pbi)
%h1 = histogram(Z_saved);
%h2 = histogram(poisson_saved);
%h = histogram2(Z_saved,beta_saved,100, 'displaystyle', 'Tile');

x = Z_saved;
y = beta_saved;
z = exp(log_posterior_pbi);
numbins=50;

[N,xedges, yedges] = histcounts2(x,y,numbins);

for i=1:numbins
    for j=1:numbins
        p = find(Z_saved>=xedges(i) & Z_saved<=xedges(i+1) & beta_saved>=yedges(j) & beta_saved<=yedges(j+1));

       if isempty(p)
           bin_max(i,j)=0;
       else
           bin_max(i,j) = mean(z(p));
       end
    end
end

%figure;
%surf(xedges(1:numbins),yedges(1:numbins),bin_max', 'EdgeColor','none', 'FaceAlpha', 0.8)
contour(xedges(1:numbins),yedges(1:numbins),bin_max', 100)
standard_figure
xlabel('Acoustic impedance')
ylabel('Vs')
colorbar