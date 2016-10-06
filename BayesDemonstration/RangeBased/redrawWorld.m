function redrawWorld(belief)
        %set(hWorldFig,'CData',belief);
        colormap(cool)
        h = bar3(-50:50,belief,0.9);
        xdat=get(h,'XData');
        lightangle(-45,30)
        for i=1:length(xdat)
            xdat{i}=xdat{i}-51*ones(size(xdat{i}));
            set(h(i),'XData',xdat{i});
        end
        set(gca,'XTick',-50:50:50);
        set(gca,'YDir','normal');
        colorbar
        for ks = 1:length(h)
            zdata = get(h(ks),'ZData');
            set(h(ks),'CData', zdata);
            set(h(ks),'FaceColor', 'interp');
        end
        axis tight
        
    end