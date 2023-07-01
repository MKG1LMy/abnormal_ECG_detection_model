%% 预处理
num_abnormal = numel(abnormal_data(:,1));
num_normal = numel(normal_data(:,1));

if num_normal/num_abnormal>3/2
F = randperm(num_normal,num_abnormal);%平衡数据规模
normal_data = normal_data(F,:);
num_normal = numel(normal_data(:,1));
end

if num_abnormal/num_normal>3/2
F = randperm(num_abnormal,num_normal);%平衡数据规模
abnormal_data = abnormal_data(F,:);
num_abnormal = numel(abnormal_data(:,1));
end

X_data = zeros(num_abnormal,3);
N_data = zeros(num_normal,3);


%载入样本
for i =1:num_abnormal
     X_data(i,1) = abnormal_data{i,3};
%      X_data(i,2) = abnormal_data{i,4};
%      X_data(i,3) = abnormal_data{i,5};
end

for i =1:num_normal
     N_data(i,1) = normal_data{i,3};
%      N_data(i,2) = normal_data{i,4};
%      N_data(i,3) = normal_data{i,5};
end



%% 训练
X_abnormal = X_data; % 异常心电图信号矩阵，每行代表一个样本
X_normal = N_data; % 正常心电图信号矩阵，每行代表一个样本

% 创建标签向量
Y_abnormal = ones(size(X_abnormal, 1), 1); % 异常心电图信号的标签为1
Y_normal = -ones(size(X_normal, 1), 1); % 正常心电图信号的标签为-1

% 合并数据和标签
X = [X_abnormal; X_normal];
Y = [Y_abnormal; Y_normal];

% 划分训练集和测试集
cv = cvpartition(size(X, 1), 'HoldOut', 0.2); % 将数据集划分为训练集和测试集，这里使用了80%的数据作为训练集
X_train = X(training(cv), :); % 训练集特征
Y_train = Y(training(cv), :); % 训练集标签
X_test = X(test(cv), :); % 测试集特征
Y_test = Y(test(cv), :); % 测试集标签

% SVM模型训练
 svmModel = fitcsvm(X_train, Y_train, 'KernelFunction', 'linear', 'Standardize', true);
% svmModel = fitclinear(X_train, Y_train);

% 预测
Y_pred = predict(svmModel, X_test);

% 计算准确率
accuracy = sum(Y_pred == Y_test) / numel(Y_test);

% 显示结果
disp(['准确率：', num2str(accuracy)]);

%% 绘图
figure;
scatter3(X_test(Y_test == 1, 1), X_test(Y_test == 1, 2), X_test(Y_test == 1, 3), 'r', 'filled'); % 异常心电图信号点为红色
hold on;
scatter3(X_test(Y_test == -1, 1), X_test(Y_test == -1, 2), X_test(Y_test == -1, 3), 'b', 'filled'); % 正常心电图信号点为蓝色
xlabel('Feature 1');
ylabel('Feature 2');
zlabel('Feature 3');
title('3D Scatter Plot of ECG Signals');
legend('Abnormal', 'Normal');
hold off;



%% 绘图2

% 预测结果可视化
figure;
hold on;
plot(Y_test, 'bo', 'MarkerFaceColor', 'b'); % 真实标签
plot(Y_pred, 'rx', 'MarkerSize', 10); % 预测标签
xlabel('Sample Index');
ylabel('Label');
legend('True Label', 'Predicted Label');
title('Prediction Results');
hold off;

%% 绘图3
% 混淆矩阵可视化
figure;
cm = confusionchart(Y_test, Y_pred);
cm.Title = 'Confusion Matrix';
cm.XLabel = 'Predicted Label';
cm.YLabel = 'True Label';

%% huitu

figure;
hold on;
gscatter(X_test(:, 1), X_test(:, 2), Y_pred, 'rb', 'xo'); % 根据预测结果绘制不同的标记和颜色
xlabel('Feature 1');
ylabel('Feature 2');
title('Decision Boundary of SVM Classifier');
legend('Abnormal', 'Normal');
hold off;

%% 12
figure;
subplot(1,3,1);
histogram(X_data(:,1), 'FaceColor', 'r'); % 异常样本的RR间期相邻最大差
hold on;
histogram(N_data(:,1), 'FaceColor', 'b'); % 正常样本的RR间期相邻最大差
xlabel('RR间期相邻最大差');
ylabel('频数');
legend('Abnormal', 'Normal');

subplot(1,3,2);
histogram(X_data(:,2), 'FaceColor', 'r'); % 异常样本的时间段内RR间期的标准差
hold on;
histogram(N_data(:,2), 'FaceColor', 'b'); % 正常样本的时间段内RR间期的标准差
xlabel('时间段内RR间期的标准差');
ylabel('频数');
legend('Abnormal', 'Normal');

subplot(1,3,3);
histogram(X_data(:,3), 'FaceColor', 'r'); % 异常样本的时间段内高频功率
hold on;
histogram(N_data(:,3), 'FaceColor', 'b'); % 正常样本的时间段内高频功率
xlabel('时间段内高频功率');
ylabel('频数');
legend('Abnormal', 'Normal');

suptitle('样本特征分布');
