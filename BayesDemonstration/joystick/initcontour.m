function hCline = initcontour(s_number,belief)
hCline = cell(s_number,1);
for ii=1:s_number
    hCline{ii,1} = contour(belief{1,ii},5);
end
end