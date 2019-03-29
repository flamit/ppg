function ROI = construct_ROI(allboxPoints,k,blank,ROItrim,ROIregion)
%constructs region of interest based on paramaters chosen
%full face, forehead, forehead and cheeks

%extract box coordinates of full ROI
box_coord = allboxPoints(:,2*k-1:2*k);
a1 = box_coord(1,1); a2 = box_coord(1,2); b1 = box_coord(2,1);
b2 = box_coord(2,2); c1 = box_coord(3,1); c2 = box_coord(3,2);
d1 = box_coord(4,1); d2 = box_coord(4,2);

%construct trimmed ROI
a11 = a1+ROItrim*(b1-a1); a21 = a2+ROItrim*(b2-a2);
b11 = a1+(1-ROItrim)*(b1-a1); b21 = a2+(1-ROItrim)*(b2-a2);
c11 = c1+ROItrim*(d1-c1); c21 = c2+ROItrim*(d2-c2);
d11 = c1+(1-ROItrim)*(d1-c1); d21 = c2+(1-ROItrim)*(d2-c2);

trimmed = insertShape(blank,'FilledPolygon',...
        [a11,a21,b11,b21,c11,c21,d11,d21],...
        'Color', 'white', 'Opacity', 1);
    
%construct eyeless and mouthless ROI
mouthless = insertShape(blank,'FilledPolygon',...
            [a11,a21,b11,b21,...
            b11+0.7*(c11-b11),b21+0.7*(c21-b21),...
            a11+0.7*(d11-a11),a21+0.7*(d21-a21)],...
            'Color', 'white', 'Opacity', 1);
        
eyelessmouthless = get_shape(mouthless,a11,a21,b11,b21,d11,d21,...
    0,0.2,0.4,0.5,'black');
eyelessmouthless = get_shape(eyelessmouthless,a11,a21,b11,b21,d11,d21,...
    0.6,0.2,1,0.5,'black');

%construct small forehead ROI
forehead = get_shape(blank,a11,a21,b11,b21,d11,d21,...
    0.3,0,0.7,0.2,'white');
forehead = get_shape(forehead,a11,a21,b11,b21,d11,d21,...
    0.4,0.2,0.6,0.3,'white');

%construct larger forehead and cheeks
forecheek = get_shape(blank,a11,a21,b11,b21,d11,d21,...
    0.2,0,0.8,0.2,'white');
forecheek = get_shape(forecheek,a11,a21,b11,b21,d11,d21,...
    0.4,0.2,0.6,0.3,'white');
forecheek = get_shape(forecheek,a11,a21,b11,b21,d11,d21,...
    0.1,0.5,0.3,0.7,'white');
forecheek = get_shape(forecheek,a11,a21,b11,b21,d11,d21,...
    0.7,0.5,0.9,0.7,'white');

    %trimmed face ROI
if ROIregion == 1
    ROI = trimmed;
    %eyeless and mouthless ROI
elseif ROIregion == 2
    ROI = eyelessmouthless;
    %forehead ROI
elseif ROIregion == 3
    ROI = forecheek;
    %forehead and cheek ROI
elseif ROIregion == 4
    ROI = forehead;
    
end
end