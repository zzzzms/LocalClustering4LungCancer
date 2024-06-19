% ====================================================================== %
% Local Cluster Extraction (LCE) Algorihtm on Lung Cancer dataset.
% ====================================================================== %

% ========================= Acknowledgement =============================
% This code is written by Mr. Ye Tian under Dr. Ming-Jun Lai's supervision 
% in March 2024. It is further modified by Dr. Ming-Jun Lai and Dr. Zhaiming Shen
% in June 2024.
% =======================================================================


clear, close all, clc, warning off

load("weightedLungCancerAdj_LN2.mat");
uAw = LN_Aw;
Ind = 1:1097; %Ind = Ind(rand);
Ind_benign = find(Ind >= 1 & Ind <= 120);
Ind_malignant = find(Ind >= 121 & Ind <= 681);
Ind_normal = find(Ind >= 682 & Ind <= 1097);

N=100; %number of repetitions

score_benign_benign_mat = zeros(N,1);
score_benign_malignant_mat = zeros(N,1);
score_benign_normal_mat = zeros(N,1);
score_malignant_benign_mat = zeros(N,1);
score_malignant_malignant_mat = zeros(N,1);
score_malignant_normal_mat = zeros(N,1);
score_normal_benign_mat = zeros(N,1);
score_normal_malignant_mat = zeros(N,1);
score_normal_normal_mat = zeros(N,1);
score_mat = zeros(N,1);

tic
for k=1:N
    
    benign_sample = datasample(Ind_benign, 10, 'Replace', false);
    malignant_sample = datasample(Ind_malignant, 50, 'Replace', false);
    normal_sample = datasample(Ind_normal, 40, 'Replace', false);
    sample = [benign_sample,malignant_sample,normal_sample];
    
    for i = 1:100
        A = uAw;
        tempsample = sample;
        Gamma = tempsample(i);
        n0 = 6; epsilon = 0.8;

        Cluster = main_CS_LCE(A,Gamma,n0,epsilon,3,0.0);
        GammaIndex = find(Cluster == Gamma, 1);
        Cluster(GammaIndex) = [];

        n1 = sum(Cluster < 121);
        n2 = sum(Cluster >= 121 & Cluster <= 681);
        n3 = sum(Cluster > 681);
        
        if i<=10
            if n1 >= max(n2,n3)
                score_benign_benign_mat(k) = score_benign_benign_mat(k) + 1;
            elseif n2 >= max(n1,n3)
                score_benign_malignant_mat(k) = score_benign_malignant_mat(k) + 1;
            elseif n3 >= max(n1,n2)
                score_benign_normal_mat(k) = score_benign_normal_mat(k) + 1;
            end
        elseif (i > 10 && i<= 60)
            if n1 >= max(n2,n3)
                score_malignant_benign_mat(k) = score_malignant_benign_mat(k) + 1;
            elseif n2 >= max(n1,n3)
                score_malignant_malignant_mat(k) = score_malignant_malignant_mat(k) + 1;
            elseif n3 >= max(n1,n2)
                score_malignant_normal_mat(k) = score_malignant_normal_mat(k) + 1;
            end
        elseif (i > 60)
            if n1 >= max(n2,n3)
                score_normal_benign_mat(k) = score_normal_benign_mat(k) + 1;
            elseif n2 >= max(n1,n3)
                score_normal_malignant_mat(k) = score_normal_malignant_mat(k) + 1;
            elseif n3 >= max(n1,n2)
                score_normal_normal_mat(k) = score_normal_normal_mat(k) + 1;
            end
        end
    end
end
toc

score_benign = mean(score_benign_benign_mat)/.1
score_malignant = mean(score_malignant_malignant_mat)/.5
score_normal = mean(score_normal_normal_mat)/.4
score = mean(score_benign_benign_mat) + mean(score_malignant_malignant_mat) + mean(score_normal_normal_mat);

STD = [std(score_benign_benign_mat/0.1), std(score_malignant_malignant_mat/0.5), std(score_normal_normal_mat/0.4), std(score_benign_benign_mat+score_malignant_malignant_mat+score_normal_normal_mat)]

disp(['The accuracy of our identification for benign images is ',num2str(score_benign),'%'])
disp(['The accuracy of our identification for malignant images is ',num2str(score_malignant),'%'])
disp(['The accuracy of our identification for normal images is ',num2str(score_normal),'%'])
disp(['The overall accuracy of our identification is ',num2str(score),'%'])
disp(['over ', num2str(N), ' repetitive tests '])

confusion_matrix = [[mean(score_benign_benign_mat), mean(score_benign_malignant_mat), mean(score_benign_normal_mat)],
    [mean(score_malignant_benign_mat), mean(score_malignant_malignant_mat), mean(score_malignant_normal_mat)],
    [mean(score_normal_benign_mat), mean(score_normal_malignant_mat), mean(score_normal_normal_mat)]]