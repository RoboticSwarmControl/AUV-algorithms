function h = drawPass(belief)
        %set(hWorldFig,'CData',belief);
        colormap(parula)
        h = surf(belief,'EdgeColor','interp');
%         camlight headlight;
%         lighting gouraud;
% %         h.FaceLighting = 'gouraud';
%          h.AmbientStrength = 0.1;
%          h.DiffuseStrength = 0.5;
%          h.SpecularStrength = 0.5;
%          h.SpecularExponent = 10;
%          h.BackFaceLighting = 'unlit';
%         colorbar
        for ks = 1:length(h)
            zdata = get(h(ks),'ZData');
            set(h(ks),'CData', zdata);
            set(h(ks),'FaceColor', 'interp');
        end
        axis tight
       
        
    end