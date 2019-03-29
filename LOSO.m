function [final_TP,final_FP,final_FN,final_TN,overall_ACC,...
    F1_1,F1_2,overall_F1,overall_MCC] = LOSO(V,P,F,All_inputs)
% leave-one-subject out cross validation, feature selection, Naive Bayes

thresh = 0.3; %Fisher LD threshold
fishers = zeros(P,F);
prediction_matrix = zeros(V,P);
overall_accuracies = zeros(5,P); %ACC,TP,FP,FN,TN

for subj = 1:P    
    %test data
    test_set_inputs = All_inputs(:,2:F+1,subj);
    test_set_labels = All_inputs(:,1,subj);
    
    %train data setup
    train_setup_inputs = All_inputs;
    train_setup_inputs(:,:,subj) = [];
    train_set_inputs = zeros(V*(P-1),F);
    train_set_labels = zeros(V*(P-1),1);
    
    %reshape training data
    for sub = 1:P-1
        train_set_inputs((sub-1)*V+1:sub*V,:) = train_setup_inputs(:,2:F+1,sub);
        train_set_labels((sub-1)*V+1:sub*V) = train_setup_inputs(:,1,sub);
    end
    
    %normalise
    train_means = mean(train_set_inputs,1);
    train_stds = std(train_set_inputs,0,1);
    
    train_set_inputs = (train_set_inputs-train_means)./train_stds;
    test_set_inputs = (test_set_inputs-train_means)./train_stds;
    
    %feature selection
    train_lows = train_set_inputs(logical(~train_set_labels),:);
    train_highs = train_set_inputs(logical(train_set_labels),:);
        
    J_Fisher = abs(mean(train_lows,1)-mean(train_highs,1))...
        ./(var(train_lows,0,1)+var(train_highs,0,1));
    
    fishers(subj,:) = J_Fisher;
    
    %select features above a threshold
    chosen_features = logical(J_Fisher>thresh);
    
    %trim training and test set inputs accordingly
    train_lows = train_lows(:,chosen_features);
    train_highs = train_highs(:,chosen_features);
    test_set_inputs = test_set_inputs(:,chosen_features);
        
    %Train Gaussian Naive Bayes
    p1 = sum(train_set_labels)/length(train_set_labels);
    p0 = 1-p1;
    
    Gauss_mean0 = mean(train_lows,1);
    Gauss_var0 = var(train_lows,0,1);
    
    Gauss_mean1 = mean(train_highs,1);
    Gauss_var1 = var(train_highs,0,1);
    
    %Inference
    comparison = zeros(V,2);
    ground_truth = test_set_labels;
    
    comparison(:,1) = log(p0)+sum(-0.5*log(2*pi*Gauss_var0)...
        -0.5*Gauss_var0.^-1.*(test_set_inputs-Gauss_mean0).^2,2);
    
    comparison(:,2) = log(p1)+sum(-0.5*log(2*pi*Gauss_var1)...
        -0.5*Gauss_var1.^-1.*(test_set_inputs-Gauss_mean1).^2,2);
    
    [~,prediction] = max(comparison,[],2);
    prediction = prediction - 1;
    
    prediction_matrix(:,subj) = prediction; 

    %evaulation
    overall_accuracies(1,subj) = sum(logical(ground_truth==prediction))/V; %ACC
    outcome = prediction + 2*ground_truth;
    TP = sum(logical(outcome == 3));
    FP = sum(logical(outcome == 1));
    FN = sum(logical(outcome == 2));
    TN = sum(logical(outcome == 0));
    
    overall_accuracies(2,subj) = TP; %TP
    overall_accuracies(3,subj) = FP; %FP
    overall_accuracies(4,subj) = FN; %FN
    overall_accuracies(5,subj) = TN; %TN
      
end

%Fisher's linear discriminant bar chart and histogram
fishers(isnan(fishers))=0;
fish_bar = mean(fishers,1);
figure
bar(fish_bar)
ylim([0,1])
hold on
plot(thresh*ones(1,F),'r-')
xlabel('feature')
ylabel('mean Fisher Linear Discriminant')
print -depsc LOSOFisher.eps

%ACC, F1-SCORE, MCC AND CONFUSION MATRIX
final_TP = sum(overall_accuracies(2,:));
final_FP = sum(overall_accuracies(3,:));
final_FN = sum(overall_accuracies(4,:));
final_TN = sum(overall_accuracies(5,:));
overall_ACC = (final_TP+final_TN)/(final_TP+final_TN+final_FP+final_FN);
%overall F1 score
if final_FP == 0
    precision = 1;
else
    precision = final_TP/(final_TP+final_FP);
end
if final_FN == 0 || final_TP == 0
    recall = 1;
else
    recall = final_TP/(final_TP+final_FN);
end
%switch around
TP2 = final_TN;
FP2 = final_FN;
FN2 = final_FP;
if FP2 == 0
    precision2 = 1;
else
    precision2 = TP2/(TP2+FP2);
end
if FN2 == 0 || TP2 == 0
    recall2 = 1;
else
    recall2 = TP2/(TP2+FN2);
end
F1_1 = 2*precision*recall/(precision+recall); %first F1-score
F1_2 = 2*precision2*recall2/(precision2+recall2); %second F1-score
overall_F1 = (F1_1+F1_2)/2; %average F1-score

%Matthews Correlation Coefficient
overall_MCC = (final_TP*final_TN - final_FP*final_FN)/...
    sqrt((final_TP+final_FP)*(final_TP+final_FN)...
    *(final_TN+final_FP)*(final_TN+final_FN));

end
