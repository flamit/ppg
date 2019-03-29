%Naive Bayes stress detection classifier
close all
clear
fileread = 'GroundTruth_features.xlsx';
F = 7; %number of features
P = 16; %number of participants

%% 1. low stress = maths easy, rest easy, rest hard
% load in the features
V = 8; %number of segments
All_inputs = zeros(V,F+1,P);

for p = 1:P
    start = 10*(p-1) + 5;
    features = xlsread(fileread,'features',...
        strcat('M',num2str(start),':T',num2str(start+V-1)));
    All_inputs(:,:,p) = features;
end

% leave-one-subject-out cross validation, feature selection, Naive Bayes
[final_TP_1,final_FP_1,final_FN_1,final_TN_1,overall_ACC_1,...
    F1_1_1,F1_2_1,overall_F1_1,overall_MCC_1] = LOSO(V,P,F,All_inputs)

%% 2. low stress = maths easy
% load in the features
V = 4; %number of segments
All_inputs = zeros(V,F+1,P);

for p = 1:P
    start = 10*(p-1) + 5;
    %maths 2
    features1 = xlsread(fileread,'features',...
        strcat('M',num2str(start),':T',num2str(start+1)));
    %maths 4
    features2 = xlsread(fileread,'features',...
        strcat('M',num2str(start+4),':T',num2str(start+5)));
    
    All_inputs(:,:,p) = [features1;features2];
end

% leave-one-subject-out cross validation, feature selection, Naive Bayes
[final_TP_2,final_FP_2,final_FN_2,final_TN_2,overall_ACC_2,...
    F1_1_2,F1_2_2,overall_F1_2,overall_MCC_2] = LOSO(V,P,F,All_inputs)

%% 3. low stress = rest hard
% load in the features
V = 4; %number of segments
All_inputs = zeros(V,F+1,P);

for p = 1:P
    start = 10*(p-1) + 5;
    pot_features = xlsread(fileread,'features',...
        strcat('M',num2str(start),':T',num2str(start+2*V-1)));
    %if hard maths task first
    if pot_features(1,1) == 1
        features = pot_features(1:4,:);
    %if hard maths task second
    else
        features = pot_features(5:8,:);
    end
    
    All_inputs(:,:,p) = features;
end

% leave-one-subject-out cross validation, feature selection, Naive Bayes
[final_TP_3,final_FP_3,final_FN_3,final_TN_3,overall_ACC_3,...
    F1_1_3,F1_2_3,overall_F1_3,overall_MCC_3] = LOSO(V,P,F,All_inputs)

%% 4. low stress = rest easy
% load in the features
V = 4; %number of segments
All_inputs = zeros(V,F+1,P);

for p = 1:P
    start = 10*(p-1) + 5;
    pot_features = xlsread(fileread,'features',...
        strcat('M',num2str(start),':T',num2str(start+2*V-1)));
    %if hard maths task first
    if pot_features(1,1) == 1
        features1 = pot_features(1:2,:);
        features2 = pot_features(7:8,:);        
    %if hard maths task second
    else
        features1 = pot_features(5:6,:);
        features2 = pot_features(3:4,:);
    end
    
    All_inputs(:,:,p) = [features1;features2];
end

% leave-one-subject-out cross validation, feature selection, Naive Bayes
[final_TP_4,final_FP_4,final_FN_4,final_TN_4,overall_ACC_4,...
    F1_1_4,F1_2_4,overall_F1_4,overall_MCC_4] = LOSO(V,P,F,All_inputs)

%% Export all to a spreadsheet
export = 0;
if export == 1
    starting_cell = 'B54';
    confusion = [final_TP_1,final_FP_1,final_FN_1,final_TN_1,overall_ACC_1,...
        F1_1_1,F1_2_1,overall_F1_1,overall_MCC_1;...
        final_TP_2,final_FP_2,final_FN_2,final_TN_2,overall_ACC_2,...
        F1_1_2,F1_2_2,overall_F1_2,overall_MCC_2;...
        final_TP_3,final_FP_3,final_FN_3,final_TN_3,overall_ACC_3,...
        F1_1_3,F1_2_3,overall_F1_3,overall_MCC_3;...
        final_TP_4,final_FP_4,final_FN_4,final_TN_4,overall_ACC_4,...
        F1_1_4,F1_2_4,overall_F1_4,overall_MCC_4];
    
    xlswrite('confusion.xlsx',confusion,1,starting_cell);
end





    
    




