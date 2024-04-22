function [h1, h2, h3, h4, h5, h6, h7, h8, fig1, fig2] = plot_hists(models_pbi, poisson_saved, Z_saved)
        fig1 = figure;
        subplot(2,3,1);
        h1=histogram(models_pbi(:,1), 'BinWidth',10);
        xlabel('{rho} 1 (kg/m3)');
        ylabel('freq')
        xlim([0 4000]);
        standard_figure
        pbaspect([1 1 1])
        
        subplot(2,3,2);
        h2 = histogram(models_pbi(:,3), 'BinWidth',10);
        xlabel('Vp 1 (m/s)');
        ylabel('freq')
        xlim([0 8000]);
        standard_figure
                pbaspect([1 1 1])


        subplot(2,3,3);
        h3=histogram(models_pbi(:,5), 'BinWidth',10);
        xlabel('Vs 1 (m/s)');
        ylabel('freq')
        xlim([0 3000]);
        standard_figure
        pbaspect([1 1 1])

        
        subplot(2,3,4)
        h4=histogram(models_pbi(:,2), 'BinWidth',25);
        xlabel('{rho} 2 (kg/m3)');
        ylabel('freq')
        xlim([0 4000]);
        standard_figure
                pbaspect([1 1 1])

        subplot(2,3,5)
        h5=histogram(models_pbi(:,4), 'BinWidth',25);
        xlabel('Vp 2 (m/s)');
        ylabel('freq')
        xlim([0 8000]);
        standard_figure
        
        pbaspect([1 1 1])


        subplot(2,3,6)
        h6=histogram(models_pbi(:,6), 'BinWidth',25);
        xlabel('Vs 2 (m/s)');
        ylabel('freq')
        xlim([0 3000]);
        standard_figure
        pbaspect([1 1 1])

        set(gcf,'Name','histograms_density-Vp-Vs');
                set(gcf, 'position', [400, 400, 1400, 800])


        fig2 = figure;         
        set(gcf, 'position', [400, 400, 1400, 800])

        subplot(1,2,1)

        if min(poisson_saved) < 0 && max(poisson_saved) > 0.5
            edges = [min(poisson_saved), 0:0.0001:0.5, max(poisson_saved)];
        else
            edge_grid = 0:0.005:0.5;
            min_edge = find(edge_grid<=min(poisson_saved), 1, 'last');
            max_edge = find(edge_grid>=max(poisson_saved), 1, 'first');
            edges = edge_grid(min_edge:max_edge);

        end
        h7=histogram(poisson_saved,edges);
        xlim([0,0.5]);
        xlabel('Poisson ratio');
        ylabel('freq');
        standard_figure
        pbaspect([1 1 1])

        subplot(1,2,2)
        h8=histogram(Z_saved, 'BinWidth',100000);
        xlabel('Acoustic impedance (kg m^{-2} s^{-1})');
        ylabel('freq');
        standard_figure
        pbaspect([1 1 1])

        set(gcf,'Name','histograms_Z-poisson');

        
        

