%% 读取数据
data_num = [100:109,111:119,121:124,200:203,205,207:210,212:215,217,219:223,228,230:234];
records = cell(1, 48);  % 创建一个大小为48的单元格数组来存储数据
annotations = cell(1, 48);  % 创建一个大小为48的单元格数组来存储注释信息

for i = 1:48
    recordname = num2str(data_num(i));  % 将编号转换为字符串形式
    records{i} = rdsamp(recordname);  % 读取数据并存储在对应的单元格中
    annotations{i} = rdann(recordname, 'atr');  % 读取注释信息并存储在对应的单元格中
end

for i = 1:48
    records{i} = records{i}(:, 1:2);  % 删除最后一列
end

sampling_rate = 360;  % 采样率


%% HRV
for i = 1:48
    record_data = records{i}(:, 2);  % 获取当前记录的ECG数据

    % R峰检测
    [~, R_peaks_indices] = findpeaks(record_data, 'MinPeakHeight', mean(record_data), 'MinPeakDistance', sampling_rate/2);

    % 计算每两个R峰之间的采样点距离间隔
    RR_intervals = diff(R_peaks_indices);
    
    %计算R峰数量及R峰间隔数量
    num_RR_intervals = numel(RR_intervals);
    num_R_peaks_indices = numel(R_peaks_indices);

    %创建待储存HRV
    RR_records = zeros(size(records{i}(:,1)));

    %填入数据

    RR_records(R_peaks_indices(1)) = RR_intervals(1);

    for a =1:num_RR_intervals
        RR_records(R_peaks_indices(a+1)) = RR_intervals(a) ;
    end
    
    
    records{i}(:,3) = RR_records;

%     填补起始和终止区段
%     RR_records(1:R_peaks_indices(1)) = RR_intervals(1);
%     RR_records(R_peaks_indices(num_R_peaks_indices):numel(records{i}(:,1))) = RR_intervals(num_RR_intervals);
% 
%     %填补R峰间数据
%     for a = 1:num_R_peaks_indices-1
%         RR_records(R_peaks_indices(a):R_peaks_indices(a+1)) = RR_intervals(a);
%     end
% 
%     records{i}(:,3) = RR_records;
% 
end





%% 分割数据
time_length=30;
segment_length = time_length*sampling_rate;  % 每个片段的长度

segments = struct('data', cell(48, 1),'ann',cell(48,1),'hrv',cell(48,1),'sc',cell(48,1),'sc2',cell(48,1),'sc3',cell(48,1));  % 创建新的结构体数组用于存储分割后的数据和标识位

for i = 1:48
    num_segments = floor(size(records{i}, 1) / segment_length);  % 计算每组数据可以分割成的片段数
    
    segments(i).data = cell(1, num_segments);  % 创建单元格数组用于存储每个片段的数据
    
    for j = 1:num_segments
        start_index = (j-1) * segment_length + 1;  % 起始索引
        end_index = j * segment_length;  % 结束索引
        
        segments(i).data{j}(:,2) = records{i}(start_index:end_index, 2);  % 存储每个片段的数据
        segments(i).data{j}(:,1) = records{i}(start_index:end_index, 1);  % 存储每个片段的标识位
        segments(i).hrv{j}(:,1) = records{i}(start_index:end_index, 3); %储存HRV数据

        nzos = nonzeros(records{i}(start_index:end_index, 3));
        c_nzos = [];
        for k =1:numel(nzos)-1
            c_nzos=[c_nzos,abs(nzos(k)-nzos(k+1))];
        end
        segments(i).sc{j} = min(nzos);
        segments(i).sc2{j}=std(nzos);
        segments(i).sc3{j}=calculate_high_freq_power(records{i}(start_index:end_index, 2),sampling_rate);
    end
end

%% 去基线

% for i =1:48
%     num_segments = numel(segments(i).data);
% 
%     for j =1:num_segments
%         segments(i).data{j}(:,2) = segments(i).data{j}(:,2)-mean(segments(i).data{j}(:,2));
%     end
% end

%% 归一化

% % for i =1:48
% %     num_segments = numel(segments(i).data);
% % 
% %     for j =1:num_segments
% %         segments(i).data{j}(:,2) = segments(i).data{j}(:,2)./mean(segments(i).data{j}(:,2));
% %     end
% % end

%% 逐点标记 
% for i = 1:48
%     
%     num_segments = numel(segments(i).data);  % 获取每组数据的片段数
%     
%     for j = 1:num_segments
%         
%         segment_data = segments(i).data{j};  % 获取当前片段的数据
% 
%         for z = 1:segment_length
%             
%             segment_num = segment_data(z, 1);  % 获取当前片段的标识位
% 
%             for k = 1:numel(annotations{i})
%                 if annotations{i}(k).sampleNumber == segment_num
%                 annotation_type(z,1) = annotations{i}(k).typeMnemonic;  % 获取标识位的类型
%                 segments(i).ann{j}=annotation_type;
%                 end
%             end
%         end
%     end
% end

%% 逐段标记
for i = 1:48
    num_segments = numel(segments(i).data);  % 获取每组数据的片段数
    bj = 1;

    for j = 1:num_segments
    
        segment_data = segments(i).data{j};  % 获取当前片段的数据
        annotation_type = 'N';
        segments(i).ann{j} = 'N';
        
        for z = 1:segment_length
           
            segment_num = segment_data(z, 1);  % 获取当前片段的标识位
            
            for k = bj:numel(annotations{i})
                
                if annotations{i}(k).sampleNumber == segment_num
                    annotation_type =annotations{i}(k).typeMnemonic;  % 获取标识位的类型
                    bj = k;
                    break;
                end

                if annotations{i}(k).sampleNumber > segment_num
                    break;
                end

            end

            if annotation_type ~='N'
               segments(i).ann{j} = annotation_type;
               break;
            end

        end        
    end
end

%% 分类存储数据
% N_data = cell(48, 2);  % 存储标记为'N'的数据
% A_data = cell(48, 2);  % 存储标记为'A'的数据
% 
% for i = 1:48
%     num_segments = numel(segments(i).data);  % 获取每组数据的片段数
%     
%     N_indices = find(strcmp(segments(i).ann, 'N'));  % 获取标记为'N'的片段索引
%     A_indices = find(strcmp(segments(i).ann, 'A'));  % 获取标记为'A'的片段索引
%         
%     for j = 1:numel(N_indices)
%         N_data{i,1}(j,:) = segments(i).data{N_indices(j)}(:, 2)';  % 将标记为'N'的数据存储到相应的单元格中，每段数据占一行
%         N_data{i,2}(j,:) = segments(i).hrv{N_indices(j)}(:,1)';
%     end
%     
%     for j = 1:numel(A_indices)
%         A_data{i,1}(j,:) = segments(i).data{A_indices(j)}(:, 2)';  % 将标记为'A'的数据存储到相应的单元格中，每段数据占一行
%         A_data{i,2}(j,:) = segments(i).hrv{N_indices(j)}(:,1)';
%     end
% end

%% 分类

abnormal_data =cell(0);
normal_data=cell(0);
PIT_abnomal = 1;
PIT_normal = 1;

jl = [];
for i =1:48
    num_segments = numel(segments(i).data);
    for j = 1:num_segments
        if segments(i).ann{j} =='R'||segments(i).ann{j} =='L'||segments(i).ann{j} =='/'
            jl = [jl;i];
            break;
        end
    end
end

CC = setdiff(1:48,jl);

for i = CC
    num_segments = numel(segments(i).data);  % 获取每组数据的片段数
    
    N_indices = find(strcmp(segments(i).ann, 'N'));  % 获取标记为'N'的片段索引
    A_indices = find(strcmp(segments(i).ann, 'A')|strcmp(segments(i).ann, 'V'|strcmp(segments(i).ann, 'a')));

    for j = 1:numel(N_indices)
        normal_data{PIT_normal,1}= segments(i).data{N_indices(j)}(:, 2)';  % 将标记为'N'的数据存储到相应的单元格中，每段数据占一行
        normal_data{PIT_normal,2} = segments(i).hrv{N_indices(j)}(:,1)';
        normal_data{PIT_normal,3} = segments(i).sc{N_indices(j)};
        normal_data{PIT_normal,4} = segments(i).sc2{N_indices(j)};
        normal_data{PIT_normal,5} = segments(i).sc3{N_indices(j)};
        normal_data{PIT_normal,6} = [i,N_indices(j)];
        PIT_normal= PIT_normal+1;
    end
    
    for j = 1:numel(A_indices)
        abnormal_data{PIT_abnomal,1} = segments(i).data{A_indices(j)}(:, 2)';  % 将标记为'A'的数据存储到相应的单元格中，每段数据占一行
        abnormal_data{PIT_abnomal,2} = segments(i).hrv{A_indices(j)}(:,1)';
        abnormal_data{PIT_abnomal,3} = segments(i).sc{A_indices(j)};
        abnormal_data{PIT_abnomal,4} = segments(i).sc2{A_indices(j)};
        abnormal_data{PIT_abnomal,5} = segments(i).sc3{A_indices(j)};
        abnormal_data{PIT_abnomal,6} = [i,A_indices(j)];
        PIT_abnomal = PIT_abnomal+1;
    end
end


