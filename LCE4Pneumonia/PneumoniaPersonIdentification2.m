% ====================================================================== %
% Local Cluster Extraction (LCE) Algorihtm on Pneumonia dataset.
% ====================================================================== %

% ========================= Acknowledgement =============================
% This code is written by Mr. Ye Tian under Dr. Ming-Jun Lai's supervision 
% in March 2024. It is further modified by Dr. Ming-Jun Lai and Dr. Zhaiming Shen
% in June 2024.
% =======================================================================

clear, close all, clc, warning off
load("Pneum.mat");

% Pneumonia Images
Ind = 1:4708;
Ind_normal = find(Ind >= 1 & Ind <= 1214);
Ind_malignant = find(Ind >= 1215 & Ind <= 4708);
 
N=100; Result=[]; %N is the number of repetitive tests. 

score_malignant_malignant_mat = zeros(N,1);
score_malignant_normal_mat = zeros(N,1);
score_normal_malignant_mat = zeros(N,1);
score_normal_normal_mat = zeros(N,1);
score_mat = zeros(N,1);

tic
for k = 1:N
    k
    malignant_sample = datasample(Ind_malignant, 75, 'Replace', false);
    normal_sample = datasample(Ind_normal, 25, 'Replace', false);
    sample = [normal_sample,malignant_sample];
    
    for i = 1:100
        tempsample = sample;
        Gamma = tempsample(i);
%        tempsample(i) = [];
%        A(tempsample, :) = 0;
%        A(:, tempsample) = 0;

        n0 = 8; epsilon = 0.8;
        Cluster = main_CS_LCE(A,Gamma,n0,epsilon,3,0.0);
        GammaIndex = find(Cluster == Gamma, 1);
        Cluster(GammaIndex) = [];
        [~,Cluster_size] = size(Cluster);

        if Cluster_size == 0
            continue
            
        end

        n1 = sum(Cluster <= 1214);
        n2 = sum(Cluster >= 1215);

        if (i <= 25)
            if n1 > n2
                score_normal_normal_mat(k) = score_normal_normal_mat(k) + 1;
            else
                score_normal_malignant_mat(k) = score_normal_malignant_mat(k) + 1;
            end
        else
            if n1 < n2
                score_malignant_malignant_mat(k) = score_malignant_malignant_mat(k) + 1;   
            else
                score_malignant_normal_mat(k) = score_malignant_normal_mat(k) + 1;
            end
        end
    end
end
toc

score_normal = mean(score_normal_normal_mat)/.25
score_malignant = mean(score_malignant_malignant_mat)/.75
score = mean(score_normal_normal_mat) + mean(score_malignant_malignant_mat);

STD = [std(score_normal_normal_mat/0.25), std(score_malignant_malignant_mat/0.75), std(score_malignant_malignant_mat+score_normal_normal_mat)]

disp(['The accuracy of our identification for normal images is ',num2str(score_normal),'%'])
disp(['The accuracy of our identification for malignant images is ',num2str(score_malignant),'%'])
disp(['The overall accuracy of our identification is ',num2str(score),'%'])
disp(['over ', num2str(N), ' repetitive tests '])

confusion_matrix = [[mean(score_normal_normal_mat), mean(score_normal_malignant_mat)],
                    [mean(score_malignant_normal_mat), mean(score_malignant_malignant_mat)]]

% return 
% %Print some images of pneumonia images. 
% load Pneum.mat
% B=[];
% for j=0:9
%     A=[];
% for i=1:10
%     a=L(i+j*10+2000,:); a1=reshape(a,28,28); a2=rot90(a1,-1);
%     A=[A a2];
% end
% B=[B;A];
% end
% figure, image(B)
% colormap gray(256) 
