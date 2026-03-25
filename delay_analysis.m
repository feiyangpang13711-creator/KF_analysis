%% ========================================================================
%   1) 数据遍历：从 data_import 读取每个 c3d 记录，按文件名提取被试与条件。
%   2) 轨迹拼接：将左右手 X/Y 序列拼接为完整 300s 轨迹。
%   3) Kalman 分段拟合：对下采样后的轨迹按 segment_size 求取观测噪声 r_std_x/y。
%   4) CCG/对齐与距离：计算对齐后的逐点误差（欧氏/|X|/|Y|），
%      并给出 20s 段内平均、全程曲线；【新增】同时返回
%      有符号误差 (ex,ey) 及瞬时方向 theta=atan2(ey,ex)（弧度）。
%   5) 存储结构 results/summary：单被试层面存储分段量、均值、详细/全程序列等；
%      群体层面聚合并计算均值±SEM。
%   6) 可视化：
%       - visualizeSingleSubjectTimeseries：单被试 300s 与 20s 段内平均曲线
%       - visualizeXYDistancesAndAccumulation：群体 X/Y 20s 曲线与误差积累
%       - visualizeFull300SecondError：群体 300s 全程误差
%       - visualizeDirectionDistributionPerSubject：单被试偏移方向分布（极坐标直方图）
%
%  新增模块说明（方向分布）：
%    - 方向定义：对齐后响应相对目标的瞬时偏移向量 e = [ex, ey] = [rx - tx, ry - ty]。
%    - 方向角：theta = atan2(ey, ex)，弧度，范围 (-pi, pi]。
%    - 可视化：对每个条件绘制 polarhistogram（归一化为概率），
%              并叠加平均合成向量（宽度为结果向量范数 R，指向平均方向）。
%% ========================================================================

data = data_import;                         % 外部导入的数据结构
conditions = {'vis','low','mid','high'};
subject_ids = 101:116;

% ---- 分析参数 ----
segment_size = 20;         % Kalman 分段长度（下采样后点数）
sampling_rate = 50;        % Hz
seconds_per_segment = 20;  % detailed 段长（秒）
Q = 1*1;                   % Kalman 过程噪声
warning('off');

% 颜色与标签（与现有风格一致）
COND_KEYS = {'vis','low','mid','hig'};
COND_NAME = containers.Map({'vis','low','mid','hig'},{'视觉','近','中','远'});
C.vis = [0.2,0.6,0.4]; C.low = [0.2,0.4,0.8]; C.mid = [0.8,0.4,0.6]; C.hig = [0.9,0.6,0.3];

%% ============================ 遍历数据，生成 results ============================
results = struct();

for data_idx = 1:length(data)
    current_data = data(data_idx).c3d;
    filename     = data(data_idx).filename;
    [subject_id, condition] = extract_info_from_filename(filename);
    if ~ismember(subject_id, subject_ids) || ~ismember(condition, conditions)
        fprintf('跳过文件: %s (不符合预期格式)\n', filename{1}); 
        continue;
    end
    if strcmp(condition,'quest.kinarm')
        fprintf('跳过QUEST条件文件: %s\n', filename{1}); 
        continue;
    end

    fprintf('正在处理被试 %d, 条件 %s (文件: %s)\n', subject_id, condition, filename{1});
    subject_field = ['subject_' num2str(subject_id)];
    cond = condition(1:3);                 % vis/low/mid/hig（保持你的缩写）

    % 1) 拼接原始轨迹
    [rx,ry,lx,ly] = concat_coords(current_data, 3, 300);
    response_ori  = [rx'; ry'];
    target_ori    = [lx'; ly'];

    % 2) 分段 Kalman（等价于你原来的分段拟合）
    fprintf('执行分段卡尔曼滤波分析...\n');
    [subject_mean_x, subject_mean_y, r_std_x, r_std_y] = ...
        do_kalman_segments(target_ori, response_ori, segment_size, Q);

    % 3) CCG/对齐与距离 + 【新增】方向（有符号误差与方向角）
    fprintf('执行交叉相关轨迹分析...\n');
    [mean_distance, mean_distance_x, mean_distance_y, ...
     segment_distances, segment_distances_x, segment_distances_y, ...
     detailed_means, detailed_means_x, detailed_means_y, ...
     full_distance, full_distance_x, full_distance_y, ...
     full_ex, full_ey, full_theta] = ...
        do_ccg_align_and_distance(target_ori, response_ori, ...
                                  seconds_per_segment, sampling_rate);

    % 4) 存储（字段名与原逻辑保持一致；新增：full_ex/full_ey/error_theta）
    if ~isfield(results, subject_field), results.(subject_field) = struct(); end
    results.(subject_field).(cond).r_seg_mean_x        = subject_mean_x;
    results.(subject_field).(cond).r_seg_mean_y        = subject_mean_y;
    results.(subject_field).(cond).segments_x          = r_std_x;
    results.(subject_field).(cond).segments_y          = r_std_y;
    results.(subject_field).(cond).mean_distance       = mean_distance;
    results.(subject_field).(cond).mean_distance_x     = mean_distance_x;
    results.(subject_field).(cond).mean_distance_y     = mean_distance_y;
    results.(subject_field).(cond).segment_distances   = segment_distances;
    results.(subject_field).(cond).segment_distances_x = segment_distances_x;
    results.(subject_field).(cond).segment_distances_y = segment_distances_y;
    results.(subject_field).(cond).detailed_distances   = detailed_means;
    results.(subject_field).(cond).detailed_distances_x = detailed_means_x;
    results.(subject_field).(cond).detailed_distances_y = detailed_means_y;
    results.(subject_field).(cond).full_distance   = full_distance;
    results.(subject_field).(cond).full_distance_x = full_distance_x;
    results.(subject_field).(cond).full_distance_y = full_distance_y;
    % 新增：有符号误差与方向角
    results.(subject_field).(cond).full_ex        = full_ex;
    results.(subject_field).(cond).full_ey        = full_ey;
    results.(subject_field).(cond).error_theta    = full_theta; % 弧度

    fprintf('存储被试 %s 条件 %s 的数据完成\n', subject_field, cond);
end

%% ============================ 汇总 summary ============================
summary = struct();
for i = 1:numel(COND_KEYS)
    summary.(COND_KEYS{i}) = init_summary_bucket();
end

% 聚合
subject_fields = fieldnames(results);
for s = 1:numel(subject_fields)
    subject = subject_fields{s};
    for i = 1:numel(COND_KEYS)
        cond = COND_KEYS{i};
        if ~isfield(results.(subject), cond), continue; end
        R = results.(subject).(cond);

        summary.(cond).kalman_subjects_x = [summary.(cond).kalman_subjects_x; R.r_seg_mean_x];
        summary.(cond).kalman_subjects_y = [summary.(cond).kalman_subjects_y; R.r_seg_mean_y];

        if isfield(R,'segments_x'), summary.(cond).kalman_segments_x = [summary.(cond).kalman_segments_x; R.segments_x']; end
        if isfield(R,'segments_y'), summary.(cond).kalman_segments_y = [summary.(cond).kalman_segments_y; R.segments_y']; end

        summary.(cond).ccg_distances   = [summary.(cond).ccg_distances;   R.mean_distance];
        summary.(cond).ccg_distances_x = [summary.(cond).ccg_distances_x; R.mean_distance_x];
        summary.(cond).ccg_distances_y = [summary.(cond).ccg_distances_y; R.mean_distance_y];

        if isfield(R,'segment_distances')
            [summary.(cond).ccg_segments,   ~] = align_append(summary.(cond).ccg_segments,   R.segment_distances');
            [summary.(cond).ccg_segments_x, ~] = align_append(summary.(cond).ccg_segments_x, R.segment_distances_x');
            [summary.(cond).ccg_segments_y, ~] = align_append(summary.(cond).ccg_segments_y, R.segment_distances_y');
        end

        if isfield(R,'detailed_distances')
            summary.(cond).detailed_distances{end+1}   = R.detailed_distances;
            summary.(cond).detailed_distances_x{end+1} = R.detailed_distances_x;
            summary.(cond).detailed_distances_y{end+1} = R.detailed_distances_y;
        end
        if isfield(R,'full_distance')
            summary.(cond).full_distances{end+1}   = R.full_distance;
            summary.(cond).full_distances_x{end+1} = R.full_distance_x;
            summary.(cond).full_distances_y{end+1} = R.full_distance_y;
        end
    end
end

% 统计
for i = 1:numel(COND_KEYS)
    cond = COND_KEYS{i};
    if isempty(summary.(cond).kalman_subjects_x)
        fprintf('条件 %s 没有有效数据，跳过\n', cond); 
        continue;
    end
    S = summary.(cond);

    S.mean_kalman_x = mean(S.kalman_subjects_x); S.std_kalman_x = std(S.kalman_subjects_x);
    S.mean_kalman_y = mean(S.kalman_subjects_y); S.std_kalman_y = std(S.kalman_subjects_y);

    S.mean_ccg_distance   = mean(S.ccg_distances);   S.std_ccg_distance   = std(S.ccg_distances);
    S.mean_ccg_distance_x = mean(S.ccg_distances_x); S.std_ccg_distance_x = std(S.ccg_distances_x);
    S.mean_ccg_distance_y = mean(S.ccg_distances_y); S.std_ccg_distance_y = std(S.ccg_distances_y);

    if ~isempty(S.ccg_segments)
        S.mean_ccg_segments   = mean(S.ccg_segments,   1, 'omitnan');
        S.std_ccg_segments    = std( S.ccg_segments,   0, 1, 'omitnan');
        S.mean_ccg_segments_x = mean(S.ccg_segments_x, 1, 'omitnan');
        S.std_ccg_segments_x  = std( S.ccg_segments_x, 0, 1, 'omitnan');
        S.mean_ccg_segments_y = mean(S.ccg_segments_y, 1, 'omitnan');
        S.std_ccg_segments_y  = std( S.ccg_segments_y, 0, 1, 'omitnan');
    end

    if ~isempty(S.detailed_distances)
        n = numel(S.detailed_distances);
        pts = seconds_per_segment*sampling_rate;
        [M,MX,MY] = deal(zeros(n,pts));
        for k=1:n
            dx = S.detailed_distances{k}(:,1);  if size(dx,1)>size(dx,2), dx=dx'; end
            dy = S.detailed_distances_x{k}(:,1);if size(dy,1)>size(dy,2), dy=dy'; end
            dz = S.detailed_distances_y{k}(:,1);if size(dz,1)>size(dz,2), dz=dz'; end
            if numel(dx)==pts, M(k,:)=dx; MX(k,:)=dy; MY(k,:)=dz; end
        end
        S.avg_detailed_distances   = mean(M,  1, 'omitnan'); S.std_detailed_distances   = std(M,  0, 1, 'omitnan');
        S.avg_detailed_distances_x = mean(MX, 1, 'omitnan'); S.std_detailed_distances_x = std(MX, 0, 1, 'omitnan');
        S.avg_detailed_distances_y = mean(MY, 1, 'omitnan'); S.std_detailed_distances_y = std(MY, 0, 1, 'omitnan');
    end

    if ~isempty(S.full_distances)
        L = inf;
        for k=1:numel(S.full_distances)
            if ~isempty(S.full_distances{k}), L = min(L, numel(S.full_distances{k})); end
        end
        if isfinite(L)
            [A,AX,AY] = deal([]);
            for k=1:numel(S.full_distances)
                if ~isempty(S.full_distances{k})
                    A  = [A;  S.full_distances{k}(1:L)];
                    AX = [AX; S.full_distances_x{k}(1:L)];
                    AY = [AY; S.full_distances_y{k}(1:L)];
                end
            end
            S.avg_full_distances   = mean(A,  1, 'omitnan'); S.std_full_distances   = std(A,  0, 1, 'omitnan');
            S.avg_full_distances_x = mean(AX, 1, 'omitnan'); S.std_full_distances_x = std(AX, 0, 1, 'omitnan');
            S.avg_full_distances_y = mean(AY, 1, 'omitnan'); S.std_full_distances_y = std(AY, 0, 1, 'omitnan');
        end
    end
    summary.(cond) = S; % 回写

    fprintf('条件 %s 结果:\n', cond);
    fprintf(' 平均卡尔曼X不确定性: %.2f ± %.2f\n', S.mean_kalman_x, S.std_kalman_x);
    fprintf(' 平均卡尔曼Y不确定性: %.2f ± %.2f\n', S.mean_kalman_y, S.std_kalman_y);
    fprintf(' 平均CCG距离: %.4f ± %.4f\n', S.mean_ccg_distance, S.std_ccg_distance);
    fprintf(' 平均CCG距离X: %.4f ± %.4f\n', S.mean_ccg_distance_x, S.std_ccg_distance_x);
    fprintf(' 平均CCG距离Y: %.4f ± %.4f\n', S.mean_ccg_distance_y, S.std_ccg_distance_y);
end

%% ============================ 可视化调用 ============================
% 单被试：（300s + 20s 段内平均）
visualizeSingleSubjectTimeseries(results, subject_ids, seconds_per_segment, sampling_rate, C, COND_NAME);
%% 群体：
visualizeXYDistancesAndAccumulation(summary, COND_KEYS, seconds_per_segment, C, COND_NAME);
%% 
visualizeFull300SecondError(summary, COND_KEYS, sampling_rate, C, COND_NAME);
%% 新增：单被试偏移方向分布（极坐标直方图 + 平均向量）
visualizeDirectionDistributionPerSubject(results, subject_ids, C, COND_NAME);
%% ============================ 本地函数区 ============================

function [rx,ry,lx,ly] = concat_coords(current_data, i0, i1)
    [RX] = {current_data.Right_HandX}; [RY] = {current_data.Right_HandY};
    [LX] = {current_data.Left_HandX};  [LY] = {current_data.Left_HandY};
    rx=[]; ry=[]; lx=[]; ly=[];
    for i = i0:i1
        rx = [rx; RX{i}]; ry = [ry; RY{i}];
        lx = [lx; LX{i}]; ly = [ly; LY{i}];
    end
end

function [mx,my, r_std_x, r_std_y] = do_kalman_segments(target_ori, response_ori, segment_size, Q)
    % 下采样（每 50 点）
    target   = target_ori(:, 1:50:end);
    response = response_ori(:, 1:50:end);
    L = min(size(target,2), size(response,2));
    target = target(:,1:L); response = response(:,1:L);

    num_segments = floor(L/segment_size);
    r_std_x = zeros(num_segments,1);
    r_std_y = zeros(num_segments,1);
    opt = optimset('Display','off');
    r0  = 5;

    for seg=1:num_segments
        idx = (seg-1)*segment_size + (1:segment_size);
        tx = (target(1,idx)-mean(target(1,idx)))*100;
        ty = (target(2,idx)-mean(target(2,idx)))*100;
        rx = (response(1,idx)-mean(response(1,idx)))*100;
        ry = (response(2,idx)-mean(response(2,idx)))*100;

        % 需要外部函数 negLogLikelihoodr(r, Q, target, response)
        [rlogx,~] = fminunc(@(r) negLogLikelihoodr(r,Q,tx,rx), r0, opt);
        [rlogy,~] = fminunc(@(r) negLogLikelihoodr(r,Q,ty,ry), r0, opt);
        r_std_x(seg) = sqrt(exp(rlogx));
        r_std_y(seg) = sqrt(exp(rlogy));
    end
    % >10 过滤（与你原来一致）
    fx = r_std_x(r_std_x<=10); fy = r_std_y(r_std_y<=10);
    mx = ternary(isempty(fx), mean(r_std_x), mean(fx));
    my = ternary(isempty(fy), mean(r_std_y), mean(fy));
end

function [mean_d, mean_dx, mean_dy, ...
          seg_d, seg_dx, seg_dy, ...
          det_d, det_dx, det_dy, ...
          full_d, full_dx, full_dy, ...
          full_ex, full_ey, full_theta] = ...
          do_ccg_align_and_distance(target_ori, response_ori, seconds_per_segment, sampling_rate)
    % 功能：对齐目标与响应轨迹，计算误差（欧氏/|X|/|Y|），
    %      并返回全程逐点误差，以及 20s 段内逐点平均；
    %      【新增】返回有符号误差 (ex,ey) 与方向角 theta=atan2(ey,ex)（弧度）。

    % 转 cm
    tx = target_ori(1,:)*100; ty = target_ori(2,:)*100;
    rx = response_ori(1,:)*100; ry = response_ori(2,:)*100;
    rx = rx - 50; % 保留你原来的偏移

    % xcorr 求滞后（随后强制置 0，保持原输出）
    [xcx,lx] = xcorr(tx, rx, 'coeff'); %#ok<ASGLU>
    [xcy,ly] = xcorr(ty, ry, 'coeff'); %#ok<ASGLU>
    lag = 0; % 原逻辑：计算最优滞后后仍置 0

    if lag>0
        atx = tx(lag+1:end); arx = rx(1:end-lag);
        aty = ty(lag+1:end); ary = ry(1:end-lag);
    else
        L = abs(lag);
        atx = tx(1:end-L); arx = rx(L+1:end);
        aty = ty(1:end-L); ary = ry(L+1:end);
    end

    % 对齐后长度
    M = min([numel(atx),numel(arx),numel(aty),numel(ary)]);
    atx=atx(1:M); arx=arx(1:M); aty=aty(1:M); ary=ary(1:M);

    % 有符号误差（响应 - 目标）
    ex = (arx - atx);
    ey = (ary - aty);

    % 距离（原有的绝对/欧式）
    d  = sqrt((atx-arx).^2 + (aty-ary).^2);
    dx = abs(atx-arx); 
    dy = abs(aty-ary);

    mean_d  = mean(d);  mean_dx = mean(dx); mean_dy = mean(dy);

    % 按 20 s 段求均值
    pts = seconds_per_segment*sampling_rate;
    numSeg = floor(M/pts);
    seg_d  = zeros(numSeg,1); seg_dx = zeros(numSeg,1); seg_dy = zeros(numSeg,1);
    det_d  = zeros(1,pts);    det_dx = zeros(1,pts);    det_dy = zeros(1,pts);
    valid = 0;
    for s=1:numSeg
        J = (s-1)*pts + (1:pts);
        seg_d(s)  = mean(d(J));
        seg_dx(s) = mean(dx(J));
        seg_dy(s) = mean(dy(J));
        % 逐点累加（段内平均）
        valid = valid + 1;
        det_d  = det_d  + d(J)';
        det_dx = det_dx + dx(J)';
        det_dy = det_dy + dy(J)';
    end
    if valid>0
        det_d  = det_d / valid; det_dx = det_dx / valid; det_dy = det_dy / valid;
    end

    % 全程（用于单被试看 300s）
    full_d  = d;
    full_dx = dx;
    full_dy = dy;

    % 新增：返回有符号误差与方向角（弧度）
    full_ex    = ex;
    full_ey    = ey;
    full_theta = atan2(ey, ex); % (-pi, pi]
end

function S = init_summary_bucket()
    S = struct( ...
        'kalman_subjects_x',[], 'kalman_subjects_y',[], ...
        'kalman_segments_x',[], 'kalman_segments_y',[], ...
        'ccg_distances',[], 'ccg_distances_x',[], 'ccg_distances_y',[], ...
        'ccg_segments',[], 'ccg_segments_x',[], 'ccg_segments_y',[], ...
        'detailed_distances',{cell(0)}, 'detailed_distances_x',{cell(0)}, 'detailed_distances_y',{cell(0)}, ...
        'full_distances',{cell(0)}, 'full_distances_x',{cell(0)}, 'full_distances_y',{cell(0)} );
end

function [acc, cur_out] = align_append(acc, cur)
    % 把行向量 cur 追加到 acc（对齐 min 长度）
    cur_out = cur;
    if isempty(acc)
        acc = cur(:)'; return;
    end
    m = min(size(acc,2), numel(cur));
    if size(acc,2)>m, acc = acc(:,1:m); end
    if numel(cur)>m,  cur_out = cur(1:m);   end
    acc = [acc; cur_out(:)'];
end

function r = ternary(cond, a, b)
    if cond, r=a; else, r=b; end
end

function [subject_id, condition] = extract_info_from_filename(filename)
    subject_id = []; condition = [];
    if iscell(filename), filename = filename{1}; end
    tokens = regexp(filename, '^(\d+)_([^.]+)', 'tokens', 'once');
    if ~isempty(tokens)
        subject_id = str2double(tokens{1});
        condition  = tokens{2};  % e.g., 'vis.kinarm'
    end
end

%% ============================ 可视化：新增 单被试 ============================
function visualizeSingleSubjectTimeseries(results, subject_ids, seconds_per_segment, sampling_rate, C, COND_NAME)
    cond_fields = {'vis','low','mid','hig'};

    % 选择被试
    avail = {}; 
    for id = subject_ids
        f = ['subject_' num2str(id)]; 
        if isfield(results,f), avail{end+1}=num2str(id); end
    end
    if isempty(avail), errordlg('没有可用被试'); return; end
    [ix,ok] = listdlg('PromptString','选择一个被试','ListString',avail,'SelectionMode','single');
    if ~ok, return; end
    sid = str2double(avail{ix}); sf = ['subject_' num2str(sid)];

    % 该被试可用条件
    conds = {};
    for i=1:numel(cond_fields)
        if isfield(results.(sf),cond_fields{i}), conds{end+1}=cond_fields{i}; end
    end
    if isempty(conds), errordlg('该被试无条件数据'); return; end

    %% 图1：300s 全程（欧式/X/Y）
    figure('Position',[100,100,1200,800]); tl = tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
    title(tl, ['被试 ' num2str(sid) ' · 300s误差随时间'], 'FontWeight','bold');

    nexttile; hold on; leg1={};
    for i=1:numel(conds)
        cond = conds{i};
        if ~isfield(results.(sf).(cond),'full_distance'), continue; end
        d  = results.(sf).(cond).full_distance;
        t  = (0:numel(d)-1)/sampling_rate;
        plot(t, d, 'Color', C.(cond), 'LineWidth', 1.8); 
        leg1{end+1} = COND_NAME(cond);
    end
    grid on; xlabel('时间 (s)'); ylabel('欧式距离 (cm)'); if ~isempty(leg1), legend(leg1,'Location','best'); end; title('欧式距离');

    nexttile; hold on; leg2={};
    for i=1:numel(conds)
        cond = conds{i};
        if ~isfield(results.(sf).(cond),'full_distance_x'), continue; end
        dx = results.(sf).(cond).full_distance_x; t=(0:numel(dx)-1)/sampling_rate;
        plot(t, dx, 'Color', C.(cond), 'LineWidth', 1.8); 
        leg2{end+1} = COND_NAME(cond);
    end
    grid on; xlabel('时间 (s)'); ylabel('X距离 (cm)'); if ~isempty(leg2), legend(leg2,'Location','best'); end; title('X方向距离');

    nexttile; hold on; leg3={};
    for i=1:numel(conds)
        cond = conds{i};
        if ~isfield(results.(sf).(cond),'full_distance_y'), continue; end
        dy = results.(sf).(cond).full_distance_y; t=(0:numel(dy)-1)/sampling_rate;
        plot(t, dy, 'Color', C.(cond), 'LineWidth', 1.8); 
        leg3{end+1} = COND_NAME(cond);
    end
    grid on; xlabel('时间 (s)'); ylabel('Y距离 (cm)'); if ~isempty(leg3), legend(leg3,'Location','best'); end; title('Y方向距离');

    %% 图2：20s 段内逐点平均（欧式/X/Y；横轴 0–20s）
    pts = seconds_per_segment*sampling_rate; tx = linspace(0, seconds_per_segment, pts);
    figure('Position',[150,150,1200,800]); tl2 = tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
    title(tl2, ['被试 ' num2str(sid) ' · 20s段内逐点平均'], 'FontWeight','bold');

    nexttile; hold on; leg4={};
    for i=1:numel(conds)
        cond = conds{i};
        if ~isfield(results.(sf).(cond),'detailed_distances') || isempty(results.(sf).(cond).detailed_distances), continue; end
        dd = results.(sf).(cond).detailed_distances; if size(dd,1)>size(dd,2), dd=dd'; end
        plot(tx, dd, 'Color', C.(cond), 'LineWidth', 1.8); leg4{end+1} = COND_NAME(cond);
    end
    grid on; xlabel('时间 (s)'); ylabel('欧式距离 (cm)'); if ~isempty(leg4), legend(leg4,'Location','best'); end; title('欧式距离（段内平均）');

    nexttile; hold on; leg5={};
    for i=1:numel(conds)
        cond = conds{i};
        if ~isfield(results.(sf).(cond),'detailed_distances_x') || isempty(results.(sf).(cond).detailed_distances_x), continue; end
        ddx = results.(sf).(cond).detailed_distances_x; if size(ddx,1)>size(ddx,2), ddx=ddx'; end
        plot(tx, ddx, 'Color', C.(cond), 'LineWidth', 1.8); leg5{end+1} = COND_NAME(cond);
    end
    grid on; xlabel('时间 (s)'); ylabel('X距离 (cm)'); if ~isempty(leg5), legend(leg5,'Location','best'); end; title('X方向（段内平均）');

    nexttile; hold on; leg6={};
    for i=1:numel(conds)
        cond = conds{i};
        if ~isfield(results.(sf).(cond),'detailed_distances_y') || isempty(results.(sf).(cond).detailed_distances_y), continue; end
        ddy = results.(sf).(cond).detailed_distances_y; if size(ddy,1)>size(ddy,2), ddy=ddy'; end
        plot(tx, ddy, 'Color', C.(cond), 'LineWidth', 1.8); leg6{end+1} = COND_NAME(cond);
    end
    grid on; xlabel('时间 (s)'); ylabel('Y距离 (cm)'); if ~isempty(leg6), legend(leg6,'Location','best'); end; title('Y方向（段内平均）');
end

%% ============================ 可视化：群体（两套图） ============================
function visualizeXYDistancesAndAccumulation(summary, cond_fields, seconds_per_segment, C, COND_NAME)
    figure('Position', [100, 100, 1200, 900]);

    % 1. X方向 20s 段内均值±SEM
    subplot(2,2,1); hold on; legend_labels = {}; plot_handles = [];
    accumulation_data_x = struct('low',[],'mid',[],'hig',[]);
    for i = 1:length(cond_fields)
        cond = cond_fields{i};
        if ~isfield(summary, cond) || ~isfield(summary.(cond), 'avg_detailed_distances_x') || isempty(summary.(cond).avg_detailed_distances_x)
            continue; 
        end
        detailed_means = summary.(cond).avg_detailed_distances_x;
        detailed_stds  = summary.(cond).std_detailed_distances_x;
        num_subjects = 0;
        if isfield(summary.(cond), 'detailed_distances_x')
            subjects_data = summary.(cond).detailed_distances_x;
            if iscell(subjects_data), num_subjects = length(subjects_data); end
        end
        detailed_sems = ternary(num_subjects>0, detailed_stds/sqrt(num_subjects), detailed_stds);
        time_points = linspace(0, seconds_per_segment, length(detailed_means));
        h = plot(time_points, detailed_means, 'Color', C.(cond), 'LineWidth', 2); plot_handles=[plot_handles,h];
        fill([time_points, fliplr(time_points)], ...
             [detailed_means+detailed_sems, fliplr(detailed_means-detailed_sems)], ...
             C.(cond), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        legend_labels{end+1} = COND_NAME(cond);

        % 误差积累（low/mid/hig）
        if ismember(cond,{'low','mid','hig'}) && isfield(summary.(cond), 'detailed_distances_x') && iscell(summary.(cond).detailed_distances_x)
            samples_per_second = length(detailed_means) / seconds_per_segment;
            idx_15s = round(15*samples_per_second) + 1;
            idx_16s = round(16*samples_per_second);
            idx_18s = round(18*samples_per_second) + 1;
            idx_19s = round(19*samples_per_second);
            for subj = 1:length(summary.(cond).detailed_distances_x)
                if isempty(summary.(cond).detailed_distances_x{subj}), continue; end
                subj_data = summary.(cond).detailed_distances_x{subj};
                if size(subj_data,1) > size(subj_data,2), subj_data = subj_data'; end
                if length(subj_data) >= idx_19s
                    mean_15_16 = mean(subj_data(idx_15s:idx_16s));
                    mean_18_19 = mean(subj_data(idx_18s:idx_19s));
                    accumulation_data_x.(cond) = [accumulation_data_x.(cond), (mean_15_16 - mean_18_19)];
                end
            end
        end
    end
    xlabel('时间 (秒)'); ylabel('X方向平均距离 (cm)'); 
    if ~isempty(legend_labels) && ~isempty(plot_handles), legend(plot_handles, legend_labels, 'Location', 'best'); end
    grid on; title('X方向距离随时间变化');

    % 2. Y方向 20s 段内均值±SEM
    subplot(2,2,2); hold on; legend_labels = {}; plot_handles = [];
    accumulation_data_y = struct('low',[],'mid',[],'hig',[]);
    for i = 1:length(cond_fields)
        cond = cond_fields{i};
        if ~isfield(summary, cond) || ~isfield(summary.(cond), 'avg_detailed_distances_y') || isempty(summary.(cond).avg_detailed_distances_y)
            continue; 
        end
        detailed_means = summary.(cond).avg_detailed_distances_y;
        detailed_stds  = summary.(cond).std_detailed_distances_y;
        num_subjects = 0;
        if isfield(summary.(cond), 'detailed_distances_y')
            subjects_data = summary.(cond).detailed_distances_y;
            if iscell(subjects_data), num_subjects = length(subjects_data); end
        end
        detailed_sems = ternary(num_subjects>0, detailed_stds/sqrt(num_subjects), detailed_stds);
        time_points = linspace(0, seconds_per_segment, length(detailed_means));
        h = plot(time_points, detailed_means, 'Color', C.(cond), 'LineWidth', 2); plot_handles=[plot_handles,h];
        fill([time_points, fliplr(time_points)], ...
             [detailed_means+detailed_sems, fliplr(detailed_means-detailed_sems)], ...
             C.(cond), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        legend_labels{end+1} = COND_NAME(cond);

        % 误差积累（low/mid/hig）
        if ismember(cond,{'low','mid','hig'}) && isfield(summary.(cond), 'detailed_distances_y') && iscell(summary.(cond).detailed_distances_y)
            samples_per_second = length(detailed_means) / seconds_per_segment;
            idx_15s = round(15*samples_per_second) + 1;
            idx_16s = round(16*samples_per_second);
            idx_18s = round(18*samples_per_second) + 1;
            idx_19s = round(19*samples_per_second);
            for subj = 1:length(summary.(cond).detailed_distances_y)
                if isempty(summary.(cond).detailed_distances_y{subj}), continue; end
                subj_data = summary.(cond).detailed_distances_y{subj};
                if size(subj_data,1) > size(subj_data,2), subj_data = subj_data'; end
                if length(subj_data) >= idx_19s
                    mean_15_16 = mean(subj_data(idx_15s:idx_16s));
                    mean_18_19 = mean(subj_data(idx_18s:idx_19s));
                    accumulation_data_y.(cond) = [accumulation_data_y.(cond), (mean_15_16 - mean_18_19)];
                end
            end
        end
    end
    xlabel('时间 (秒)'); ylabel('Y方向平均距离 (cm)'); 
    if ~isempty(legend_labels) && ~isempty(plot_handles), legend(plot_handles, legend_labels, 'Location', 'best'); end
    grid on; title('Y方向距离随时间变化');

    % 3. X方向 误差积累（条形图）
    subplot(2,2,3); 
    conditions = {'近','中','远'};
    mean_accum_x = [mean(accumulation_data_x.low), mean(accumulation_data_x.mid), mean(accumulation_data_x.hig)];
    sem_accum_x  = [std(accumulation_data_x.low)/sqrt(max(1,length(accumulation_data_x.low))), ...
                    std(accumulation_data_x.mid)/sqrt(max(1,length(accumulation_data_x.mid))), ...
                    std(accumulation_data_x.hig)/sqrt(max(1,length(accumulation_data_x.hig)))];
    b = bar(1:3, mean_accum_x); hold on; b.FaceColor = 'flat';
    b.CData = [C.low; C.mid; C.hig];
    errorbar(1:3, mean_accum_x, sem_accum_x, 'k.', 'LineWidth', 1.5);
    ylabel('X方向误差积累 (cm)'); xlabel('距离条件'); xticks(1:3); xticklabels(conditions); grid on; title('X方向误差积累分析');
    % ANOVA
    accumulation_values_x = [accumulation_data_x.low, accumulation_data_x.mid, accumulation_data_x.hig];
    accumulation_groups_x = [repmat({'近'},1,length(accumulation_data_x.low)), ...
                             repmat({'中'},1,length(accumulation_data_x.mid)), ...
                             repmat({'远'},1,length(accumulation_data_x.hig))];
    if ~isempty(accumulation_values_x)
        p_x = anova1(accumulation_values_x, accumulation_groups_x, 'off');
        text(0.5, 0.9*max(eps+mean_accum_x), sprintf('ANOVA: p = %.3f', p_x), 'HorizontalAlignment', 'left');
        fprintf('\nX方向误差积累ANOVA: p = %.3f\n', p_x);
    end

    % 4. Y方向 误差积累（条形图）
    subplot(2,2,4); 
    mean_accum_y = [mean(accumulation_data_y.low), mean(accumulation_data_y.mid), mean(accumulation_data_y.hig)];
    sem_accum_y  = [std(accumulation_data_y.low)/sqrt(max(1,length(accumulation_data_y.low))), ...
                    std(accumulation_data_y.mid)/sqrt(max(1,length(accumulation_data_y.mid))), ...
                    std(accumulation_data_y.hig)/sqrt(max(1,length(accumulation_data_y.hig)))];
    b = bar(1:3, mean_accum_y); hold on; b.FaceColor = 'flat';
    b.CData = [C.low; C.mid; C.hig];
    errorbar(1:3, mean_accum_y, sem_accum_y, 'k.', 'LineWidth', 1.5);
    ylabel('Y方向误差积累 (cm)'); xlabel('距离条件'); xticks(1:3); xticklabels(conditions); grid on; title('Y方向误差积累分析');
    accumulation_values_y = [accumulation_data_y.low, accumulation_data_y.mid, accumulation_data_y.hig];
    accumulation_groups_y = [repmat({'近'},1,length(accumulation_data_y.low)), ...
                             repmat({'中'},1,length(accumulation_data_y.mid)), ...
                             repmat({'远'},1,length(accumulation_data_y.hig))];
    if ~isempty(accumulation_values_y)
        p_y = anova1(accumulation_values_y, accumulation_groups_y, 'off');
        text(0.5, 0.9*max(eps+mean_accum_y), sprintf('ANOVA: p = %.3f', p_y), 'HorizontalAlignment', 'left');
        fprintf('Y方向误差积累ANOVA: p = %.3f\n', p_y);
    end
end

function visualizeFull300SecondError(summary, cond_fields, sampling_rate, C, COND_NAME)
    figure('Position', [150, 150, 1200, 900]);

    % 1. 欧式距离
    subplot(3,1,1); hold on; legend_labels = {}; plot_handles = [];
    for i = 1:length(cond_fields)
        cond = cond_fields{i};
        if ~isfield(summary, cond) || ~isfield(summary.(cond), 'avg_full_distances') || isempty(summary.(cond).avg_full_distances)
            continue;
        end
        avg_distances = summary.(cond).avg_full_distances;
        std_distances = summary.(cond).std_full_distances;
        num_subjects = 0;
        if isfield(summary.(cond), 'full_distances')
            subjects_data = summary.(cond).full_distances;
            if iscell(subjects_data), num_subjects = length(subjects_data); end
        end
        sem_distances = ternary(num_subjects>0, std_distances/sqrt(num_subjects), std_distances);

        % downsample 渲染
        downsampling = 100;
        if length(avg_distances) > 1000
            avg_distances = avg_distances(1:downsampling:end);
            sem_distances = sem_distances(1:downsampling:end);
        end
        t = (0:length(avg_distances)-1)/sampling_rate;
        h = plot(t, avg_distances, 'Color', C.(cond), 'LineWidth', 2); plot_handles=[plot_handles,h];
        if length(avg_distances) < 1000
            fill([t, fliplr(t)], [avg_distances+sem_distances, fliplr(avg_distances-sem_distances)], ...
                C.(cond), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
        legend_labels{end+1} = COND_NAME(cond);
    end
    xlabel('时间 (秒)'); ylabel('欧式距离 (cm)'); 
    if ~isempty(legend_labels) && ~isempty(plot_handles), legend(plot_handles, legend_labels, 'Location', 'best'); end
    grid on; title('300秒内的欧式距离误差');

    % 2. X方向
    subplot(3,1,2); hold on; legend_labels = {}; plot_handles = [];
    for i = 1:length(cond_fields)
        cond = cond_fields{i};
        if ~isfield(summary, cond) || ~isfield(summary.(cond), 'avg_full_distances_x') || isempty(summary.(cond).avg_full_distances_x)
            continue;
        end
        avg_distances = summary.(cond).avg_full_distances_x;
        std_distances = summary.(cond).std_full_distances_x;
        num_subjects = 0;
        if isfield(summary.(cond), 'full_distances_x')
            subjects_data = summary.(cond).full_distances_x;
            if iscell(subjects_data), num_subjects = length(subjects_data); end
        end
        sem_distances = ternary(num_subjects>0, std_distances/sqrt(num_subjects), std_distances);
        downsampling = 100;
        if length(avg_distances) > 1000
            avg_distances = avg_distances(1:downsampling:end);
            sem_distances = sem_distances(1:downsampling:end);
        end
        t = (0:length(avg_distances)-1)/sampling_rate;
        h = plot(t, avg_distances, 'Color', C.(cond), 'LineWidth', 2); plot_handles=[plot_handles,h];
        if length(avg_distances) < 1000
            fill([t, fliplr(t)], [avg_distances+sem_distances, fliplr(avg_distances-sem_distances)], ...
                C.(cond), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
        legend_labels{end+1} = COND_NAME(cond);
    end
    xlabel('时间 (秒)'); ylabel('X方向距离 (cm)'); 
    if ~isempty(legend_labels) && ~isempty(plot_handles), legend(plot_handles, legend_labels, 'Location', 'best'); end
    grid on; title('300秒内的X方向距离误差');

    % 3. Y方向
    subplot(3,1,3); hold on; legend_labels = {}; plot_handles = [];
    for i = 1:length(cond_fields)
        cond = cond_fields{i};
        if ~isfield(summary, cond) || ~isfield(summary.(cond), 'avg_full_distances_y') || isempty(summary.(cond).avg_full_distances_y)
            continue;
        end
        avg_distances = summary.(cond).avg_full_distances_y;
        std_distances = summary.(cond).std_full_distances_y;
        num_subjects = 0;
        if isfield(summary.(cond), 'full_distances_y')
            subjects_data = summary.(cond).full_distances_y;
            if iscell(subjects_data), num_subjects = length(subjects_data); end
        end
        sem_distances = ternary(num_subjects>0, std_distances/sqrt(num_subjects), std_distances);
        downsampling = 100;
        if length(avg_distances) > 1000
            avg_distances = avg_distances(1:downsampling:end);
            sem_distances = sem_distances(1:downsampling:end);
        end
        t = (0:length(avg_distances)-1)/sampling_rate;
        h = plot(t, avg_distances, 'Color', C.(cond), 'LineWidth', 2); plot_handles=[plot_handles,h];
        if length(avg_distances) < 1000
            fill([t, fliplr(t)], [avg_distances+sem_distances, fliplr(avg_distances-sem_distances)], ...
                C.(cond), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
        legend_labels{end+1} = COND_NAME(cond);
    end
    xlabel('时间 (秒)'); ylabel('Y方向距离 (cm)'); 
    if ~isempty(legend_labels) && ~isempty(plot_handles), legend(plot_handles, legend_labels, 'Location', 'best'); end
    grid on; title('300秒内的Y方向距离误差');
end

%% ============================ 可视化：单被试方向分布 ============================
function visualizeDirectionDistributionPerSubject(results, subject_ids, C, COND_NAME)
    % 选择被试
    avail = {}; 
    for id = subject_ids
        f = ['subject_' num2str(id)]; 
        if isfield(results,f), avail{end+1}=num2str(id); end
    end
    if isempty(avail), errordlg('没有可用被试'); return; end
    [ix,ok] = listdlg('PromptString','选择一个被试(方向分布)','ListString',avail,'SelectionMode','single');
    if ~ok, return; end
    sid = str2double(avail{ix}); sf = ['subject_' num2str(sid)];

    cond_fields = {'vis','low','mid','hig'};
    conds = {};
    for i=1:numel(cond_fields)
        if isfield(results.(sf),cond_fields{i}) && isfield(results.(sf).(cond_fields{i}),'error_theta')
            conds{end+1} = cond_fields{i}; %#ok<AGROW>
        end
    end
    if isempty(conds), errordlg('该被试无方向数据'); return; end

    % 绘图：每个条件一个极坐标直方图 + 平均向量
    n = numel(conds);
    nrow = ceil(n/2); ncol = min(2,n);
    figure('Position',[200,200,1200,500+200*(nrow-1)]);
    t = tiledlayout(nrow, ncol, 'Padding','compact','TileSpacing','compact');
    title(t, ['被试 ' num2str(sid) ' · 偏移方向分布（atan2(ey,ex)）'], 'FontWeight','bold');

    nbins = 36; % 10° 分辨率
    for i=1:n
        cond = conds{i};
        nexttile; 
        theta = results.(sf).(cond).error_theta; % 弧度
        theta = theta(~isnan(theta) & isfinite(theta));
        if isempty(theta)
            text(0.5,0.5,'无有效方向数据','HorizontalAlignment','center'); axis off; continue;
        end
        % 归一化概率极坐标直方图
        polarhistogram(theta, nbins, 'Normalization','probability', 'FaceColor', C.(cond), 'FaceAlpha', 0.5);

        % 叠加平均合成向量
        Rcos = mean(cos(theta)); Rsin = mean(sin(theta));
        R = hypot(Rcos, Rsin); phi = atan2(Rsin, Rcos);
        hold on;
        % 在当前极坐标轴上画平均向量（从 0 到 R）
        polarplot([phi phi], [0 0.1], 'k-', 'LineWidth', 12*R);

        title([COND_NAME(cond) '（R=' sprintf('%.2f',R) '）']);
    end
end
