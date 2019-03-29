function roi = get_roi(allboxPoints,k,blank,ROItrim,...
    hei,wid,seg_height,seg_width)
%returns roi for selected box segment and frame (k)

%extract box coordinates of full ROI
box_coord = allboxPoints(:,2*k-1:2*k);
a1 = box_coord(1,1); a2 = box_coord(1,2); b1 = box_coord(2,1);
b2 = box_coord(2,2); c1 = box_coord(3,1); c2 = box_coord(3,2);
d1 = box_coord(4,1); d2 = box_coord(4,2);
%construct trimmed ROI
a11 = a1+ROItrim*(b1-a1); a21 = a2+ROItrim*(b2-a2);
b11 = a1+(1-ROItrim)*(b1-a1); b21 = a2+(1-ROItrim)*(b2-a2);
d11 = c1+(1-ROItrim)*(d1-c1); d21 = c2+(1-ROItrim)*(d2-c2);

block_width = (b11-a11)/seg_width;
block_width2 = (b21-a21)/seg_width;
block_height = (d11-a11)/seg_height;
block_height2 = (d21-a21)/seg_height;

seg_a1 = a11 + (wid-1)*block_width + (hei-1)*block_height;
seg_a2 = a21 + (wid-1)*block_width2 + (hei-1)*block_height2;
seg_b1 = a11 + wid*block_width + (hei-1)*block_height;
seg_b2 = a21 + wid*block_width2 + (hei-1)*block_height2;
seg_d1 = a11 + (wid-1)*block_width + hei*block_height;
seg_d2 = a21 + (wid-1)*block_width2 + hei*block_height2;
seg_c1 = a11 + wid*block_width + hei*block_height;
seg_c2 = a21 + wid*block_width2 + hei*block_height2;
roi = insertShape(blank,'FilledPolygon',...
    [seg_a1,seg_a2,seg_b1,seg_b2,seg_c1,seg_c2,seg_d1,seg_d2],...
    'Color', 'white', 'Opacity', 1);
end