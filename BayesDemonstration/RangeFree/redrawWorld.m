function redrawWorld(belief)
        %set(hWorldFig,'CData',belief);
        colormap(cool)
        h = bar3(belief,0.9);
        xdat=get(h,'XData');
        lightangle(-45,30)
        
        set(gca,'YDir','normal');
        colorbar
        for ks = 1:length(h)
            zdata = get(h(ks),'ZData');
            set(h(ks),'CData', zdata);
            set(h(ks),'FaceColor', 'interp');
        end
        axis tight
        
    end