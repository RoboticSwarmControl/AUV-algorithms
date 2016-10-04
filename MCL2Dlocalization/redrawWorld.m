function redrawWorld(belief)
        %set(hWorldFig,'CData',belief);
        colormap(jet)
        h = bar3(belief,.9);
        set(gca,'YDir','normal');
        colorbar
        for ks = 1:length(h)
            zdata = get(h(ks),'ZData');
            set(h(ks),'CData', zdata);
            set(h(ks),'FaceColor', 'interp');
        end
    end