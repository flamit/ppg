function filereads = feat_files

%list of all files to run through when extracting features
filereads = {'1_R1Time0101P1.aviROI_FS_coloursFIX',...
    '2_METime0854P1.aviROI_FS_coloursFIX',...
    '3_R3Time1528P1.aviROI_FS_coloursFIX',...
    '4_MHTime2302P1.aviROI_FS_coloursFIX',...
    '5_R5Time2938P1.aviROI_FS_coloursFIX',...
    '1_R1Time0034P2.aviROI_FS_coloursFIX',...
    '2_MHTime0901P2.aviROI_FS_coloursFIX',...
    '3_R3Time1530P2.aviROI_FS_coloursFIX',...
    '4_METime2229P2.aviROI_FS_coloursFIX',...
    '5_R5Time2831P2.aviROI_FS_coloursFIX',...
    '1_R1Time0039P3.aviROI_FS_coloursFIX',...
    '2_METime1022P3.aviROI_FS_coloursFIX',...
    '3_R3Time1620P3.aviROI_FS_coloursFIX',...
    '4_MHTime2316P3.aviROI_FS_coloursFIX',...
    '5_R5Time2954P3.aviROI_FS_coloursFIX',...
    '1_R1Time0116P4.aviROI_FS_coloursFIX',...
    '2_MHcutby1sP4.aviROI_FS_coloursFIX',...
    '3_R3Time1624P4.aviROI_FS_coloursFIX',...
    '4_METime2320P4.aviROI_FS_coloursFIX',...
    '5_R5Time3103P4.aviROI_FS_coloursFIX',...
    '1_R1Time0058P5.aviROI_FS_coloursFIX',...
    '2_METime0947P5.aviROI_FS_coloursFIX',...
    '3_R3cutby1sP5.aviROI_FS_coloursFIX',...
    '4_MHTime2248P5.aviROI_FS_coloursFIX',...
    '5_R5Time2926P5.aviROI_FS_coloursFIX',...
    '1_R1Time0154P6.aviROI_FS_coloursFIX',...
    '2_MHTime1128P6.aviROI_FS_coloursFIX',...
    '3_R3Time1809P6.aviROI_FS_coloursFIX',...
    '4_METime2607P6.aviROI_FS_coloursFIX',...
    '5_R5Time3251P6.aviROI_FS_coloursFIX',...
    '1_R1Time0103P7.aviROI_FS_coloursFIX',...
    '2_METime0746P7.aviROI_FS_coloursFIX',...
    '3_R3Time1344P7.aviROI_FS_coloursFIX',...
    '4_MHTime2030P7.aviROI_FS_coloursFIX',...
    '5_R5Time2630P7.aviROI_FS_coloursFIX',...
    '1_R1Time0358P8.aviROI_FS_coloursFIX',...
    '2_MHTime1207P8.aviROI_FS_coloursFIX',...
    '3_R3Time1855P8.aviROI_FS_coloursFIX',...
    '4_METime2700P8.aviROI_FS_coloursFIX',...
    '5_R5Time3404P8.aviROI_FS_coloursFIX',...
    '1_R1cutby1sP9.aviROI_FS_coloursFIX',...
    '2_METime1030P9.aviROI_FS_coloursFIX',...
    '3_R3Time1645P9.aviROI_FS_coloursFIX',...
    '4_MHTime2320P9.aviROI_FS_coloursFIX',...
    '5_R5Time2917P9.aviROI_FS_coloursFIX',...
    '1_R1Time0049P10.aviROI_FS_coloursFIX',...
    '2_MHTime0734P10.aviROI_FS_coloursFIX',...
    '3_R3Time1326P10.aviROI_FS_coloursFIX',...
    '4_METime2121P10.aviROI_FS_coloursFIX',...
    '5_R5Time2739P10.aviROI_FS_coloursFIX',...
    '1_R1Time0112P14_CUTVERSION.aviROI_FS_coloursFIX',...
    '2_MHTime0935P14.aviROI_FS_coloursFIX',...
    '3_R3Time1624P14.aviROI_FS_coloursFIX',...
    '4_METime2400P14.aviROI_FS_coloursFIX',...
    '5_R5Time3020P14.aviROI_FS_coloursFIX',...
    '1_R1Time0140P15.aviROI_FS_coloursFIX',...
    '2_METime0815P15.aviROI_FS_coloursFIX',...
    '3_R3Time1523P15.aviROI_FS_coloursFIX',...
    '4_MHTime2358P15.aviROI_FS_coloursFIX',...
    '5_R5Time3049P15.aviROI_FS_coloursFIX',...
    '1_R1Time0055P16.aviROI_FS_coloursFIX',...
    '2_MHTime1010P16.aviROI_FS_coloursFIX',...
    '3_R3Time1647P16.aviROI_FS_coloursFIX',...
    '4_METime2420P16.aviROI_FS_coloursFIX',...
    '5_R5Time3019P16.aviROI_FS_coloursFIX',...
    '1_R1Time0117P17.aviROI_FS_coloursFIX',...
    '2_METime0905P17.aviROI_FS_coloursFIX',...
    '3_R3Time1556P17.aviROI_FS_coloursFIX',...
    '4_MHTime2333P17.aviROI_FS_coloursFIX',...
    '5_R5Time3009P17.aviROI_FS_coloursFIX',...
    '1_R1Time0121P18.aviROI_FS_coloursFIX',...
    '2_MHTime0824P18.aviROI_FS_coloursFIX',...
    '3_R3Time1545P18.aviROI_FS_coloursFIX',...
    '4_METime2240P18.aviROI_FS_coloursFIX',...
    '5_R5Time2836P18.aviROI_FS_coloursFIX',...
    '1_R1Time0047P19.aviROI_FS_coloursFIX',...
    '2_METime0725P19.aviROI_FS_coloursFIX',...
    '3_R3Time1337P19.aviROI_FS_coloursFIX',...
    '4_MHTime2030P19.aviROI_FS_coloursFIX',...
    '5_R5Time2630P19.aviROI_FS_coloursFIX',...
    '1_R1Time0129P20.aviROI_FS_coloursFIX',...
    '2_MHTime0900P20.aviROI_FS_coloursFIX',...
    '3_R3Time1508P20.aviROI_FS_coloursFIX',...
    '4_METime2148P20.aviROI_FS_coloursFIX',...
    '5_R5Time2740P20.aviROI_FS_coloursFIX',...
    '1_R1Time0040P21.aviROI_FS_coloursFIX',...
    '2_METime0800P21.aviROI_FS_coloursFIX',...
    '3_R3Time1443P21.aviROI_FS_coloursFIX',...
    '4_MHTime2229P21.aviROI_FS_coloursFIX',...
    '5_R5Time2858P21.aviROI_FS_coloursFIX',...
    '1_R1Time0035P22.aviROI_FS_coloursFIX',...
    '2_MHTime0833P22.aviROI_FS_coloursFIX',...
    '3_R3Time1445P22.aviROI_FS_coloursFIX',...
    '4_METime2130P22.aviROI_FS_coloursFIX',...
    '5_R5Time2808P22.aviROI_FS_coloursFIX',...
    '1_R1Time0135P23.aviROI_FS_coloursFIX',...
    '2_METime0931P23.aviROI_FS_coloursFIX',...
    '3_R3Time1652P23.aviROI_FS_coloursFIX',...
    '4_MHTime2415P23.aviROI_FS_coloursFIX',...
    '5_R5Time3154P23.aviROI_FS_coloursFIX'};

end