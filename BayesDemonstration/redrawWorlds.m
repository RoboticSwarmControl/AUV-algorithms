function redrawWorlds(belief)
        %set(hWorldFig,'CData',belief);
        colormap(cool)
        h = surf(belief,'EdgeColor','none');
        camlight headlight;
        lighting gouraud;
%         h.FaceLighting = 'gouraud';
         h.AmbientStrength = 0.1;
         h.DiffuseStrength = 0.5;
         h.SpecularStrength = 0.5;
         h.SpecularExponent = 10;
         h.BackFaceLighting = 'unlit';
        colorbar
        for ks = 1:length(h)
            zdata = get(h(ks),'ZData');
            set(h(ks),'CData', zdata);
            set(h(ks),'FaceColor', 'interp');
        end
        axis tight
        xtick = [1 51 101];
        xticklabels = {'-50','0','50'};
        set(gca,'XTick', xtick);
        set(gca,'XTickLabel', xticklabels);
        ytick = [1 51 101];
        yticklabels = {'-50','0','50'};
        set(gca,'YTick', ytick);
        set(gca,'YTickLabel', yticklabels);
        
    end