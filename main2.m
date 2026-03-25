%% =========================
%  main2.m
%  脚本目标
%  1 逐个读取数据文件并识别被试编号与实验条件
%  2 对轨迹数据进行预处理与对齐
%  3 计算速度互相关曲线 CCG 并提取最大相关对应的时间滞后
%  4 通过极大似然估计得到卡尔曼观测噪声标准差 作为位置不确定性指标
%  5 对 QUEST 条件进行阈值后验更新与心理物理曲线拟合
%  6 输出单被试与群体层面的可视化结果
%% =========================
clearvars; clc;

%% ---------- 全局绘图风格 ----------
STYLE = viz_style("science");
set_global_style(STYLE);

%% ---------- 加载数据 ----------
% exam_load 返回结构体数组
% 每个元素对应一个数据文件 通常包含
% filename 文件名
% c3d      解析后的试次数据与运动学变量
data = exam_load();
data_import = data;                 % 保留一份用于后续轨迹可视化
data = data_import;

%% ---------- 固定参数 ----------
% subject_ids 被试编号范围
% conditions  实验条件列表
% feedback_indices 用于标记反馈相关的索引位置
subject_ids  = 101:116;
conditions   = {'vis.kinarm','low.kinarm','mid.kinarm','high.kinarm','quest.kinarm'};
feedback_indices = [17,18,37,38,57,58,77,78,97,98,117,118,137,138,157,158,177,178,197,198,217,218,237,238,257,258,277,278];

Fs = 50;                            % 采样率 Hz
maxLagSamp = 50;                    % 互相关最大滞后 样本点数
results = struct();                 % 汇总所有被试与条件的计算结果

%% =========================
%  主循环
%  每个文件完成一次完整分析
%  文件名解析 被试编号 条件识别
%  QUEST 条件走阈值后验与拟合流程
%  其余条件走轨迹拼接 CCG 计算 延迟对齐 下采样 卡尔曼噪声估计
%% =========================
for data_idx = 1:numel(data)
    current_data = data(data_idx).c3d;
    filename = data(data_idx).filename;
    if iscell(filename), fname = filename{1}; else, fname = filename; end

    % 文件名解析
    % 默认格式 101_vis.kinarm_xxx
    % 解析得到 subject_id 与 condition
    [subject_id, condition] = extract_info_from_filename(fname);

    % 仅处理目标被试与预设条件
    if ~ismember(subject_id, subject_ids) || ~ismember(condition, conditions)
        fprintf('跳过文件: %s (不符合预期格式)\n', fname);
        continue;
    end
    fprintf('正在处理被试 %d, 条件 %s (文件: %s)\n', subject_id, condition, fname);

    % 为每个被试建立独立字段
    % results.subject_101.vis 之类的层级结构
    subject_field = ['subject_' num2str(subject_id)];
    if ~isfield(results, subject_field), results.(subject_field) = struct(); end

    % 条件短名
    % 用于作为结构体字段名 避免出现点号导致字段名非法
    short = cond2short(condition);

    %% -------- QUEST 条件 ----------
    % QUEST 使用试次内的目标位置与被试选择更新阈值后验
    % 并对正确率随刺激强度的关系进行拟合
    if strcmp(condition,'quest.kinarm')
        quest_result = process_quest_data(current_data, subject_id, results);
        results.(subject_field).(short) = quest_result;
        continue
    end

    %% -------- 轨迹拼接 ----------
    % 将多个试次的轨迹按时间顺序竖直拼接为一段连续序列
    % Right_Hand 作为 response
    % Left_Hand  作为 target
    RX = {current_data.Right_HandX};  RY = {current_data.Right_HandY};
    LX = {current_data.Left_HandX};   LY = {current_data.Left_HandY};

    rx = []; ry = []; lx = []; ly = [];
    for i = 3:300
        rx = [rx; RX{i}];  ry = [ry; RY{i}];
        lx = [lx; LX{i}];  ly = [ly; LY{i}];
    end

    % 轨迹组织为 2 行 T 列
    % 第一行 X 第二行 Y
    response_ori = [rx'; ry'];
    target_ori   = [lx'; ly'];

    %% -------- CCG 速度互相关 ----------
    % 用一阶差分近似速度
    % 对速度进行 zscore 标准化
    % 通过 xcorr 计算互相关曲线并进行系数归一化
    targVel = zscore(diff(target_ori'))';
    respVel = zscore(diff(response_ori'))';

    ccg_x = xcorr(respVel(1,:), targVel(1,:), maxLagSamp, 'coeff');
    ccg_y = xcorr(respVel(2,:), targVel(2,:), maxLagSamp, 'coeff');

    % 取互相关最大值对应的滞后索引
    [~, idx_x] = max(ccg_x);
    [~, idx_y] = max(ccg_y);
    maxLag_x = (idx_x - (maxLagSamp+1));   % 样本点滞后 范围 -L 到 L
    maxLag_y = (idx_y - (maxLagSamp+1));

    %% -------- 延迟对齐与下采样 ----------
    % 先将两条轨迹都以起点为零点
    % 然后按 X 与 Y 维度分别进行滞后对齐
    % 最后按 Fs 步长抽取 每 50 点取 1 点 形成低频序列
    target_zero   = target_ori   - target_ori(:,1);
    response_zero = response_ori - response_ori(:,1);

    [aligned_target_x, aligned_resp_x] = lag_align(target_zero(1,:), response_zero(1,:), maxLag_x);
    [aligned_target_y, aligned_resp_y] = lag_align(target_zero(2,:), response_zero(2,:), maxLag_y);

    sampled_target_x = aligned_target_x(1:Fs:end);
    sampled_resp_x   = aligned_resp_x(1:Fs:end);
    sampled_target_y = aligned_target_y(1:Fs:end);
    sampled_resp_y   = aligned_resp_y(1:Fs:end);

    % 反馈索引
    % valid_idx 是落在当前样本长度范围内的索引
    valid_idx = feedback_indices(feedback_indices <= numel(sampled_resp_x));
    sampled_resp_x_fb = sampled_resp_x; 
    sampled_resp_y_fb = sampled_resp_y; 

    %% -------- Kalman 观测噪声估计 ----------
    % 目标序列与响应序列先去均值并裁掉开头片段
    % 以极大似然方式估计观测噪声参数 r 的对数值
    % 输出 r 的标准差 sqrt r 作为位置不确定性指标
    Q = 1;      % 过程噪声强度
    clip = 1;   % 丢弃开头的样本段
    opt = optimoptions('fminunc','Display','off');

    r_std_x = estimate_r_std(sampled_target_x, sampled_resp_x, Q, clip, opt);
    r_std_y = estimate_r_std(sampled_target_y, sampled_resp_y, Q, clip, opt);

    %% -------- 结果存储 ----------
    % 每个条件保存
    % ccg_x ccg_y 互相关曲线
    % max_lag_x max_lag_y 最大相关对应的滞后 秒
    % r_std_x r_std_y   观测噪声标准差 cm
    results.(subject_field).(short) = struct( ...
        'ccg_x', ccg_x, ...
        'ccg_y', ccg_y, ...
        'max_lag_x', maxLag_x/Fs, ...
        'max_lag_y', maxLag_y/Fs, ...
        'r_std_x', r_std_x, ...
        'r_std_y', r_std_y ...
    );
end

%% =========================
%  单被试轨迹图
%  图一 1D 随时间变化
%  每列对应一个条件 上排是 X 下排是 Y
%  target 与 response 使用不同颜色便于比较
%
%  图二 2D 轨迹
%  2 行 2 列布局 每个子图一个条件
%  response 轨迹颜色随时间由浅到深 起点与终点做标记
%% =========================
plot_subject_trajectories(data_import, subject_ids, STYLE);

%% =========================
%  群体可视化
%  CCG 群体均值曲线并显示 SEM 阴影
%  不确定性 使用散点 误差线与均值条形
%  QUEST 展示阈值估计与拟合参数
%% =========================
plot_group_ccg(results, STYLE);
plot_group_uncertainty(results, STYLE);
plot_group_quest(results, STYLE);

%% =========================
%  单被试汇总分析
%  同时展示该被试的 CCG 曲线 不确定性条形 以及 QUEST 后验分布图
%% =========================
analyze_single_subject(results, STYLE);

%% ======================================================================
%                               函数区
%% ======================================================================

function short = cond2short(condition)
% 条件名转为短字段名
% 视觉 vis
% 距离近 low
% 距离中 mid
% 距离远 hig
% QUEST que
    switch condition
        case 'vis.kinarm',   short = 'vis';
        case 'low.kinarm',   short = 'low';
        case 'mid.kinarm',   short = 'mid';
        case 'high.kinarm',  short = 'hig';
        case 'quest.kinarm', short = 'que';
        otherwise,           short = 'unk';
    end
end

function STYLE = viz_style(preset)
% 绘图风格配置
% 用于集中管理字体 字号 线宽 坐标轴线宽 网格透明度
% C 为调色板 行对应不同颜色
% ccg_cols unc_cols 分别用于互相关与不确定性图的默认配色
    if nargin<1, preset="nature"; end
    STYLE.font = "Arial";
    STYLE.font_size = 11;
    STYLE.lw = 1.6;
    STYLE.ax_lw = 1.0;
    STYLE.grid_alpha = 0.12;

    switch lower(string(preset))
        case "science"
            STYLE.C = [0.20 0.40 0.70; 0.85 0.37 0.01; 0.18 0.63 0.17; 0.55 0.35 0.64; 0.35 0.35 0.35];
        otherwise
            STYLE.C = [0.12 0.47 0.71; 0.20 0.63 0.17; 0.80 0.33 0.10; 0.60 0.31 0.64; 0.35 0.35 0.35];
    end
    STYLE.ccg_cols = STYLE.C(1:4,:);
    STYLE.unc_cols = STYLE.C(1:4,:);
end

function set_global_style(STYLE)
% 设置全局默认绘图属性
% 影响后续所有新建 figure 与 axes
% 目的是避免每次绘图都重复设置字体与坐标轴细节
    set(groot,"defaultFigureColor","w");
    set(groot,"defaultAxesFontName",STYLE.font,"defaultAxesFontSize",STYLE.font_size);
    set(groot,"defaultLineLineWidth",STYLE.lw);
    set(groot,"defaultAxesLineWidth",STYLE.ax_lw);
    set(groot,"defaultAxesTickDir","out");
    set(groot,"defaultAxesBox","off");
end

function apply_axes_style(ax, STYLE)
% 坐标轴外观统一
% 启用网格并降低透明度
% 刻度朝外
% 去掉上右边框
    grid(ax,"on"); ax.GridAlpha = STYLE.grid_alpha;
    ax.TickDir = "out"; ax.Box = "off";
end

function [subject_id, condition] = extract_info_from_filename(filename)
% 从文件名中提取信息
% 期望文件名由下划线分隔
% 第一段是被试编号
% 第二段是条件名
    subject_id = []; condition = [];
    parts = regexp(filename,'\_','split');
    if numel(parts) < 2, return; end
    subject_id = str2double(parts{1});
    condition = parts{2};
end

function [t_aligned, r_aligned] = lag_align(t, r, lag)
% 根据滞后对齐两条序列
% lag 大于 0 表示 response 相对 target 更晚出现
% lag 小于 0 表示 response 相对 target 更早出现
% 对齐后两条序列长度一致 用于后续下采样与估计
    if lag > 0
        t_aligned = t(1:end-lag);
        r_aligned = r(lag+1:end);
    elseif lag < 0
        t_aligned = t(-lag+1:end);
        r_aligned = r(1:end+lag);
    else
        t_aligned = t;
        r_aligned = r;
    end
end

function r_std = estimate_r_std(target, response, Q, clip, opt)
% 观测噪声标准差估计
% 输入
% target   对齐并下采样后的目标序列
% response 对齐并下采样后的响应序列
% Q        过程噪声强度
% clip     剪去序列前段避免启动段影响
%
% 输出
% r_std 为观测噪声参数的标准差 单位 cm
    targetc = target(:) - mean(target(1+clip:end), 'omitnan');
    responsec = response(:) - mean(response(1+clip:end), 'omitnan');

    targetc = targetc(1+clip:end) * 100;          % 单位换算到 cm
    responsec = responsec(1+clip:end) * 100;

    r0 = 3;
    r_log = fminunc(@negLogLikelihoodr, r0, opt, Q, targetc, responsec);
    r_std = sqrt(exp(r_log));
end

%% ========================= 群体 CCG =========================
function plot_group_ccg(results, STYLE)
% 群体速度互相关图
% 每个子图对应一个条件
% 横轴 Lag 秒
% 纵轴 相关系数
% 两条线分别为 X 维与 Y 维的群体均值曲线
% 阴影为均值的标准误差 SEM
    conds  = {'vis','low','mid','hig'};
    titles = {'视觉','近','中','远'};
    lag_axis = (-50:50)/50;

    f = figure("Name","Group CCG","Units","normalized","Position",[0.05 0.55 0.52 0.38]);
    t = tiledlayout(f,2,2,"TileSpacing","compact","Padding","compact");
    title(t,"群体速度互相关 CCG",'FontWeight','bold');

    for i = 1:4
        ax = nexttile(t,i); hold(ax,"on"); apply_axes_style(ax, STYLE);
        title(ax, titles{i});
        xlabel(ax,"Lag s"); ylabel(ax,"Corr");
        xline(ax,0,"-","Color",[0 0 0 0.18]);

        [ccgx, ccgy] = collect_ccg_curves(results, conds{i});
        if isempty(ccgx)
            text(ax,0.5,0.5,"无数据","Units","normalized","HorizontalAlignment","center");
            continue
        end

        mx = mean(ccgx,1,'omitnan'); sx = std(ccgx,[],1,'omitnan')/sqrt(size(ccgx,1));
        my = mean(ccgy,1,'omitnan'); sy = std(ccgy,[],1,'omitnan')/sqrt(size(ccgy,1));

        % 阴影表示 SEM 不进入图例
        fill(ax, [lag_axis fliplr(lag_axis)], [mx-sx fliplr(mx+sx)], STYLE.C(1,:), ...
            'FaceAlpha',0.16, 'EdgeColor','none','HandleVisibility','off');
        fill(ax, [lag_axis fliplr(lag_axis)], [my-sy fliplr(my+sy)], STYLE.C(4,:), ...
            'FaceAlpha',0.16, 'EdgeColor','none','HandleVisibility','off');

        % 均值曲线
        hx = plot(ax, lag_axis, mx, "-", "Color",STYLE.C(1,:), "LineWidth",STYLE.lw*0.85);
        hy = plot(ax, lag_axis, my, "-", "Color",STYLE.C(4,:), "LineWidth",STYLE.lw*0.85);

        legend(ax, [hx hy], {"X","Y"}, "Location","best", "Box","off");
        xlim(ax,[-0.5 0.5]); ylim(ax,[-0.2 1]);
    end
end

function [ccgx, ccgy] = collect_ccg_curves(results, cond)
% 从 results 中收集某一条件下所有被试的互相关曲线
% 输出矩阵
% 每行一个被试
% 每列对应一个 lag 位置
    ccgx=[]; ccgy=[];
    fields = fieldnames(results);
    for f = 1:numel(fields)
        sf = fields{f};
        if ~isfield(results.(sf), cond), continue; end
        x = results.(sf).(cond).ccg_x;
        y = results.(sf).(cond).ccg_y;
        if isempty(x) || isempty(y), continue; end
        ccgx(end+1,:) = x(:)';
        ccgy(end+1,:) = y(:)';
    end
end

%% ========================= 群体 不确定性 =========================
function plot_group_uncertainty(results, STYLE)
% 群体位置不确定性图
% 每个子图对应一个条件
% 横轴为维度类别 X 与 Y
% 纵轴为不确定性指标 sigma cm
% 可视化元素
% 半透明条形展示均值
% 散点展示个体
% 误差线展示个体标准差
    conds  = {'vis','low','mid','hig'};
    titles = {'视觉','近','中','远'};

    f = figure("Name","Group uncertainty","Units","normalized","Position",[0.05 0.10 0.52 0.36]);
    t = tiledlayout(f,2,2,"TileSpacing","compact","Padding","compact");
    title(t,"群体位置不确定性 Kalman 观测噪声 sigma",'FontWeight','bold');

    for i = 1:4
        ax = nexttile(t,i); hold(ax,"on"); apply_axes_style(ax, STYLE);
        [sx, sy] = collect_uncertainty(results, conds{i});

        title(ax, titles{i});
        ylabel(ax,"sigma cm");
        xlim(ax,[0.5 2.5]); xticks(ax,[1 2]); xticklabels(ax,{'X','Y'});

        if isempty(sx) && isempty(sy)
            text(ax,0.5,0.5,"无数据","Units","normalized","HorizontalAlignment","center");
            continue
        end

        c0 = STYLE.unc_cols(i,:);

        % 均值条形
        if ~isempty(sx), bar(ax,1,mean(sx,'omitnan'),0.55,'FaceColor',c0,'FaceAlpha',0.22,'EdgeColor','none'); end
        if ~isempty(sy), bar(ax,2,mean(sy,'omitnan'),0.55,'FaceColor',c0*0.85,'FaceAlpha',0.22,'EdgeColor','none'); end

        % 个体散点与误差线
        if ~isempty(sx)
            j=(rand(size(sx))-0.5)*0.10;
            scatter(ax,1+j,sx,26,c0,'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
            errorbar(ax,1,mean(sx,'omitnan'),std(sx,'omitnan'), 'k','CapSize',10,'LineWidth',1.1);
        end
        if ~isempty(sy)
            j=(rand(size(sy))-0.5)*0.10;
            scatter(ax,2+j,sy,26,c0*0.85,'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
            errorbar(ax,2,mean(sy,'omitnan'),std(sy,'omitnan'), 'k','CapSize',10,'LineWidth',1.1);
        end
    end
end

function [sx, sy] = collect_uncertainty(results, cond)
% 从 results 中收集某条件下所有被试的 r_std_x 与 r_std_y
% sx sy 为列向量 每个元素对应一个被试
    sx=[]; sy=[];
    fields = fieldnames(results);
    for f = 1:numel(fields)
        sf = fields{f};
        if ~isfield(results.(sf), cond), continue; end
        sx(end+1,1) = results.(sf).(cond).r_std_x;
        sy(end+1,1) = results.(sf).(cond).r_std_y;
    end
end

%% ========================= 群体 QUEST =========================
function plot_group_quest(results, STYLE)
% 群体 QUEST 阈值汇总图
% 每个被试绘制多个指标
% posterior 阈值后验均值与中位数
% 以及三种心理物理模型拟合得到的 alpha 参数
    fields = fieldnames(results);
    subs=[]; thr_mean=[]; thr_med=[]; a_log=[]; a_wei=[]; a_nor=[];

    for f = 1:numel(fields)
        sf = fields{f};
        sid = str2double(erase(sf,'subject_'));
        if ~isfield(results.(sf),'que'), continue; end
        q = results.(sf).que;

        subs(end+1,1)=sid;
        thr_mean(end+1,1)=q.threshold_mean;
        thr_med(end+1,1)=q.threshold_median;

        a_log(end+1,1)=getfield_def(q,'fit_logistic','alpha');
        a_wei(end+1,1)=getfield_def(q,'fit_weibull','alpha'); 
        a_nor(end+1,1)=getfield_def(q,'fit_normal','alpha');  
    end

    f = figure("Name","Group QUEST thresholds","Units","normalized","Position",[0.60 0.20 0.36 0.70]);
    ax = axes(f); hold(ax,'on'); apply_axes_style(ax, STYLE);
    title(ax,"群体 QUEST 阈值估计",'FontWeight','bold');
    xlabel(ax,"Subject"); ylabel(ax,"Threshold cm");

    if isempty(subs)
        text(ax,0.5,0.5,"无QUEST数据","Units","normalized","HorizontalAlignment","center");
        return
    end

    x0=(1:numel(subs))';
    xticks(ax,x0); xticklabels(ax,string(subs)); xtickangle(ax,45);

    scatter(ax,x0-0.20,thr_mean,36,STYLE.C(2,:),'filled','MarkerEdgeColor','w');
    scatter(ax,x0-0.05,thr_med, 36,STYLE.C(4,:),'filled','MarkerEdgeColor','w');
    scatter(ax,x0+0.10,a_log,   36,STYLE.C(1,:),'filled','MarkerEdgeColor','w');
    scatter(ax,x0+0.25,a_wei,   36,STYLE.C(3,:),'filled','MarkerEdgeColor','w');
    scatter(ax,x0+0.40,a_nor,   36,STYLE.C(5,:),'filled','MarkerEdgeColor','w');

    legend(ax,{'后验均值','后验中位数','Logistic拟合','Weibull拟合','Normal拟合'}, ...
        'Location','best','Box','off');
end

function v = getfield_def(q, s1, s2)
% 安全读取结构体字段
% 若字段不存在则返回 NaN
    v = NaN;
    if isfield(q,s1) && isstruct(q.(s1)) && isfield(q.(s1),s2)
        v = q.(s1).(s2);
    end
end

function [xb,pb]=bin_mean_curve(x,y,nb)
% 将数据按 x 分箱并计算每箱的均值
% xb 为每个箱的 x 均值
% pb 为每个箱的 y 均值
    x=x(:); y=y(:);
    m=~isnan(x)&~isnan(y); x=x(m); y=y(m);
    if isempty(x), xb=[]; pb=[]; return; end
    edges=linspace(min(x),max(x),nb+1);
    xb=nan(nb,1); pb=nan(nb,1);
    for i=1:nb
        idx=x>=edges(i)&x<edges(i+1);
        if any(idx), xb(i)=mean(x(idx)); pb(i)=mean(y(idx)); end
    end
    keep=~isnan(xb)&~isnan(pb); xb=xb(keep); pb=pb(keep);
end

%% ========================= 单被试轨迹图 =========================
function plot_subject_trajectories(data_import, subject_ids, STYLE)
% 单被试轨迹展示
% 交互输入被试编号
%
% 图一 1D
% 每列一个条件 上排 X 下排 Y
% target 为深灰 response 为主色
%
% 图二 2D
% 2 行 2 列布局 每个子图一个条件
% response 轨迹用渐变增强时间方向 起点空心圆 终点实心圆
    answer = inputdlg({'输入被试编号:'}, '被试轨迹图', [1 40], {num2str(subject_ids(1))});
    if isempty(answer), return; end
    subject_id = str2double(answer{1});
    if ~ismember(subject_id, subject_ids)
        errordlg(['被试 ' num2str(subject_id) ' 不存在'], '错误'); return;
    end

    conds  = {'vis.kinarm','low.kinarm','mid.kinarm','high.kinarm'};
    cname  = {'视觉','近','中','远'};
    targCol = [0.25 0.25 0.25];        % target 深灰

    % 1D X Y 随时间
    f1 = figure('Name',['Subject ' num2str(subject_id) ' (X/Y)'], ...
        'Units','normalized','Position',[0.05 0.08 0.90 0.36], 'Color','w');
    t1 = tiledlayout(f1,2,4,'TileSpacing','compact','Padding','compact');
    title(t1, ['被试 ' num2str(subject_id) '：1D 轨迹 X Y 随时间'], 'FontWeight','bold');

    % 2D 轨迹
    f2 = figure('Name',['Subject ' num2str(subject_id) ' (2D)'], ...
        'Units','normalized','Position',[0.05 0.48 0.62 0.44], 'Color','w');
    t2 = tiledlayout(f2,2,2,'TileSpacing','compact','Padding','compact');
    title(t2, ['被试 ' num2str(subject_id) '：2D 轨迹 颜色随时间渐显'], 'FontWeight','bold');

    for k = 1:4
        [ok, rx, ry, lx, ly, time_vec] = fetch_subject_condition_xy(data_import, subject_id, conds{k});

        axX = nexttile(t1,k);   hold(axX,"on"); apply_axes_style(axX, STYLE);
        axY = nexttile(t1,k+4); hold(axY,"on"); apply_axes_style(axY, STYLE);
        ax2 = nexttile(t2,k);   hold(ax2,"on"); apply_axes_style(ax2, STYLE); axis(ax2,'equal');

        title(axX,[cname{k} ' - X']); title(axY,[cname{k} ' - Y']); title(ax2,[cname{k} ' - 2D']);

        if ~ok
            text(axX,0.5,0.5,'无数据','Units','normalized','HorizontalAlignment','center'); axis(axX,'off');
            text(axY,0.5,0.5,'无数据','Units','normalized','HorizontalAlignment','center'); axis(axY,'off');
            text(ax2,0.5,0.5,'无数据','Units','normalized','HorizontalAlignment','center'); axis(ax2,'off');
            continue
        end

        lw1 = max(0.9, STYLE.lw*0.65);      % 1D 线宽
        lw2 = max(1.0, STYLE.lw*0.70);      % 2D 线宽

        % 1D target response 对比
        plot(axX, time_vec, lx, '-', 'Color',[targCol 0.70], 'LineWidth',lw1);
        plot(axX, time_vec, rx, '-', 'Color',[STYLE.C(1,:) 0.90], 'LineWidth',lw1);
        xlabel(axX,'Time s'); ylabel(axX,'X cm');

        plot(axY, time_vec, ly, '-', 'Color',[targCol 0.70], 'LineWidth',lw1);
        plot(axY, time_vec, ry, '-', 'Color',[STYLE.C(1,:) 0.90], 'LineWidth',lw1);
        xlabel(axY,'Time s'); ylabel(axY,'Y cm');

        % 2D target response 轨迹
        hT = plot(ax2, lx, ly, '-', 'Color',[targCol 0.92], 'LineWidth',lw2);
        hR = plot(ax2, NaN, NaN, '-', 'Color',[STYLE.C(1,:) 0.95], 'LineWidth',lw2); % 仅用于图例
        [hS, hE] = plot_traj_gradient(ax2, rx, ry, STYLE.C(1,:), lw2);

        xlabel(ax2,'X cm'); ylabel(ax2,'Y cm');
        legend(ax2, [hT hR hS hE], {'Target','Response','Start','End'}, 'Location','best', 'Box','off');
    end
end

function [ok, rx, ry, lx, ly, time_vec] = fetch_subject_condition_xy(data_import, subject_id, condition)
% 从 data_import 中找到指定被试与条件对应的文件
% 输出拼接后的 response 与 target 的 X Y 以及时间轴
% 数据换算为 cm
    ok=false; rx=[]; ry=[]; lx=[]; ly=[]; time_vec=[];
    for data_idx=1:numel(data_import)
        filename=data_import(data_idx).filename;
        if iscell(filename), fname=filename{1}; else, fname=filename; end
        [sid, cond]=extract_info_from_filename(fname);
        if sid~=subject_id || ~strcmp(cond,condition), continue; end

        d=data_import(data_idx).c3d;
        RX={d.Right_HandX}; RY={d.Right_HandY};
        LX={d.Left_HandX};  LY={d.Left_HandY};

        rx=[]; ry=[]; lx=[]; ly=[];
        for i=3:300
            rx=[rx; RX{i}]; ry=[ry; RY{i}];
            lx=[lx; LX{i}]; ly=[ly; LY{i}];
        end

        rx = (rx-0.5)*100;  ry=ry*100; lx=lx*100; ly=ly*100;
        time_vec=(0:numel(rx)-1)/50;
        ok=true; return
    end
end

function [hS, hE] = plot_traj_gradient(ax, x, y, baseColor, lw)
% 2D response 轨迹渐变绘制
% 线段颜色由浅到深增强时间方向
% hS 起点空心圆
% hE 终点实心圆
    n = numel(x);
    if n < 2, hS=gobjects(1); hE=gobjects(1); return; end
    if nargin < 5 || isempty(lw), lw = 1.2; end

    nseg = min(240, n-1);
    idx  = round(linspace(1, n, nseg+1));

    light = baseColor + (1-baseColor)*0.65;   % 起始更浅

    for k = 1:nseg
        i1 = idx(k); i2 = idx(k+1);
        a  = (k-1)/(nseg-1);                  % 0 到 1
        col = (1-a)*light + a*baseColor;
        alp = 0.25 + 0.70*a;                  % 透明度随时间增强
        plot(ax, x(i1:i2), y(i1:i2), '-', 'Color',[col alp], 'LineWidth',lw, 'HandleVisibility','off');
    end

    hS = scatter(ax, x(1),  y(1),  42, 'MarkerEdgeColor',baseColor, 'MarkerFaceColor','w', 'LineWidth',1.2);
    hE = scatter(ax, x(end),y(end),42, 'MarkerEdgeColor','w',        'MarkerFaceColor',baseColor, 'LineWidth',1.0);
end

%% ========================= 单被试分析 =========================
function analyze_single_subject(results, STYLE)
% 单被试分析汇总图
% 子图 1 到 4
% 分别对应四个条件的 CCG 曲线
% 曲线上用散点标记最大相关位置
%
% 子图 5
% 展示四个条件下 X 与 Y 的位置不确定性条形图
%
% 子图 6
% 展示 QUEST 阈值后验分布图
% 并在图中给出三种拟合模型的 alpha beta 参数
    answer = inputdlg({'输入被试编号 (例如: 101):'}, '单个被试分析', [1 40], {'101'});
    if isempty(answer), return; end

    subject_id = str2double(answer{1});
    subject_field = ['subject_' num2str(subject_id)];
    if ~isfield(results, subject_field)
        errordlg(['被试 ' num2str(subject_id) ' 不存在'], '错误'); return
    end

    conds  = {'vis','low','mid','hig'};
    cname  = {'视觉','近','中','远'};
    lag_axis = (-50:50)/50;

    f = figure('Name',['Single subject ' num2str(subject_id)], 'Units','normalized', ...
        'Position',[0.12 0.10 0.76 0.78], 'Color','w');
    t = tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
    title(t, ['被试 ' num2str(subject_id) '：单被试分析汇总'], 'FontWeight','bold');

    r_std_x = nan(1,4); r_std_y = nan(1,4);

    % 1 到 4 条件 CCG
    for c = 1:4
        ax = nexttile(t,c); hold(ax,'on'); apply_axes_style(ax, STYLE);
        xlabel(ax,'Lag s'); ylabel(ax,'Corr'); title(ax,[cname{c} '：CCG']);
        xline(ax,0,'-','Color',[0 0 0 0.18]);

        if ~isfield(results.(subject_field), conds{c})
            text(ax,0.5,0.5,'无数据','Units','normalized','HorizontalAlignment','center');
            continue
        end

        ccg_x = results.(subject_field).(conds{c}).ccg_x;
        ccg_y = results.(subject_field).(conds{c}).ccg_y;

        hx = plot(ax, lag_axis, ccg_x, '-', 'Color',STYLE.C(1,:), 'LineWidth',STYLE.lw*0.85);
        hy = plot(ax, lag_axis, ccg_y, '-', 'Color',STYLE.C(4,:), 'LineWidth',STYLE.lw*0.85);

        [mx, ix] = max(ccg_x); [my, iy] = max(ccg_y);
        scatter(ax, lag_axis(ix), mx, 36, STYLE.C(1,:), 'filled', 'MarkerEdgeColor','w');
        scatter(ax, lag_axis(iy), my, 36, STYLE.C(4,:), 'filled', 'MarkerEdgeColor','w');

        legend(ax,[hx hy],{'X','Y'},'Location','best','Box','off');
        xlim(ax,[-0.5 0.5]); ylim(ax,[-0.2 1]);

        r_std_x(c) = results.(subject_field).(conds{c}).r_std_x;
        r_std_y(c) = results.(subject_field).(conds{c}).r_std_y;
    end

    % 5 不确定性
    ax5 = nexttile(t,5); hold(ax5,'on'); apply_axes_style(ax5, STYLE);
    title(ax5,'位置不确定性 sigma 各条件'); ylabel(ax5,'sigma cm');
    xlim(ax5,[0.5 4.5]); xticks(ax5,1:4); xticklabels(ax5,cname);

    if all(isnan(r_std_x)) && all(isnan(r_std_y))
        text(ax5,0.5,0.5,'没有Kalman滤波数据可用','Units','normalized','HorizontalAlignment','center');
    else
        b = bar(ax5, [r_std_x(:) r_std_y(:)], 0.80, 'FaceAlpha',0.70, 'EdgeColor','none');
        b(1).FaceColor = STYLE.C(1,:);   % X
        b(2).FaceColor = STYLE.C(4,:);   % Y
        legend(ax5,{'X','Y'},'Location','best','Box','off');
    end

    % 6 QUEST
    ax6 = nexttile(t,6); hold(ax6,'on'); apply_axes_style(ax6, STYLE);
    title(ax6,'QUEST阈值'); xlabel(ax6,'Threshold cm');
    if isfield(results.(subject_field), 'que')
        q = results.(subject_field).que;
        axes(ax6); %#ok<LAXES>
        calculate_threshold(q.quest_results, q.threshold_val, true);
        box(ax6,'off');

        tstr = {
            ['Logistic:  ', num2str(q.fit_logistic.alpha,'%.2f'), ' cm (beta=', num2str(q.fit_logistic.beta,'%.2f'), ')'];
            ['Weibull:   ', num2str(q.fit_weibull.alpha,'%.2f'),  ' cm (beta=', num2str(q.fit_weibull.beta,'%.2f'),  ')'];
            ['Normal:    ', num2str(q.fit_normal.alpha,'%.2f'),   ' cm (beta=', num2str(q.fit_normal.beta,'%.2f'),   ')']
        };
        text(ax6,0.02,0.02,tstr,'Units','normalized','VerticalAlignment','bottom', ...
            'BackgroundColor',[1 1 1 0.80],'EdgeColor',[0 0 0 0.20]);
    else
        text(ax6,0.5,0.5,'没有QUEST数据可用','Units','normalized','HorizontalAlignment','center');
        axis(ax6,'off');
    end
end

%% ========================= QUEST 处理 =========================
function quest_result = process_quest_data(data_id, subject_id, results) %#ok<INUSD>
% QUEST 数据处理流程
% 将试次按轮次组织 每轮包含多个关键试次
% 第四个试次提供刺激强度
% 最后一个试次提供被试选择与正确性
%
% 通过贝叶斯更新得到阈值后验分布
% 最终输出
% threshold_mean 阈值后验均值
% threshold_median 阈值后验中位数
% fit_logistic fit_weibull fit_normal 三类心理物理模型参数
    num_trials = numel(data_id);
    threshold_val = (1:20)/2;
    origin_X = -0.25;

    threshold_prob = normpdf(threshold_val, 4, 6);
    threshold_prob = threshold_prob / sum(threshold_prob);

    num_rounds = (num_trials - 2) / 6;
    quest_results = struct('Round', [], 'Stim_Predicted', [], 'Stim_Actual', [],...
                           'Threshold_Prob', [], 'Subject_Choice', [], 'Correct', []);

    for k = 1:num_rounds
        start_trial  = 3 + (k-1)*6;
        fourth_trial = start_trial+3;
        last_trial   = start_trial+5;

        target_X = data_id(fourth_trial).Left_HandX(end);
        stim_actual = round(abs(target_X - origin_X)*200)/2;

        if stim_actual < 0 || stim_actual > 10
            warning('第%d轮刺激强度异常: %dcm', k, stim_actual);
            continue;
        end

        try
            choice = data_id(last_trial).EVENTS.LABELS{1,2};
        catch
            warning('第%d轮末尾试次%d缺少EVENTS数据', k, last_trial);
            continue;
        end

        correct_dir = 'Right Choose';
        if (target_X - origin_X) < 0, correct_dir = 'Left Choose'; end
        is_correct = strcmp(choice, correct_dir);

        x = stim_actual;
        if is_correct
            likelihood = arrayfun(@(t) Transfer(x, t), threshold_val);
        else
            likelihood = arrayfun(@(t) 1-Transfer(x, t), threshold_val);
        end

        % 依据当前后验的中位数位置给出下一轮刺激预测值
        cumu  = cumsum(threshold_prob);
        max_idx = find(cumu >= 0.5, 1, 'first');
        stim_predicted = threshold_val(max_idx);

        % 后验更新与归一化
        threshold_prob = threshold_prob .* likelihood;
        threshold_prob = threshold_prob / sum(threshold_prob);

        quest_results(k).Round = k;
        quest_results(k).Stim_Predicted = stim_predicted;
        quest_results(k).Stim_Actual = stim_actual;
        quest_results(k).Threshold_Prob = threshold_prob;
        quest_results(k).Subject_Choice = choice;
        quest_results(k).Correct = is_correct;
    end

    [threshold_mean, threshold_median] = calculate_threshold(quest_results, threshold_val, false);

    [stim_all, resp_all] = quest_to_xy(quest_results);
    fit_logistic = fit_psychometric(stim_all, resp_all, 'logistic');
    fit_weibull  = fit_psychometric(stim_all, resp_all, 'weibull');
    fit_normal   = fit_psychometric(stim_all, resp_all, 'normal');

    quest_result = struct(...
        'quest_results', quest_results, ...
        'threshold_val', threshold_val, ...
        'threshold_mean', threshold_mean, ...
        'threshold_median', threshold_median, ...
        'fit_logistic', fit_logistic, ...
        'fit_weibull',  fit_weibull, ...
        'fit_normal',   fit_normal ...
    );
end

function [threshold_mean, threshold_median] = calculate_threshold(results, threshold_val, show_plots)
% 从最后一轮的阈值后验分布提取统计量
% threshold_MAP 最大后验位置
% threshold_mean 后验均值
% threshold_median 后验中位数
% show_plots 为 true 时绘制
% 左轴 后验分布条形图并标出 MAP mean median
% 右轴 累积分布曲线并标出 0.5 水平线
    final_prob = results(end).Threshold_Prob;
    final_prob = final_prob / sum(final_prob);

    [~, map_idx] = max(final_prob);
    threshold_MAP = threshold_val(map_idx);
    threshold_mean = sum(threshold_val .* final_prob);

    cum_prob = cumsum(final_prob);
    below = find(cum_prob < 0.5, 1, 'last');
    if isempty(below)
        threshold_median = threshold_val(1);
    else
        next = min(below+1, numel(threshold_val));
        t1 = threshold_val(below); t2 = threshold_val(next);
        p1 = cum_prob(below);      p2 = cum_prob(next);
        threshold_median = t1 + (0.5-p1)/(max(p2-p1,eps))*(t2-t1);
    end

    if show_plots
        yyaxis left;
        bar(threshold_val, final_prob, 'FaceAlpha',0.55, 'FaceColor',[0.35 0.35 0.35]); hold on;
        xline(threshold_MAP,   '--', 'LineWidth',1.3, 'Color',[0.12 0.47 0.71]);
        xline(threshold_mean,  '--', 'LineWidth',1.3, 'Color',[0.20 0.63 0.17]);
        xline(threshold_median,':',  'LineWidth',1.6, 'Color',[0.60 0.31 0.64]);
        ylabel('Posterior');
        yyaxis right;
        plot(threshold_val, cum_prob, '-', 'LineWidth',1.4, 'Color',[0 0 0 0.55]);
        yline(0.5,'--','Color',[0 0 0 0.25]);
        ylabel('Cumulative');
        grid on; xlabel('Threshold cm');
        legend({'Posterior','MAP','Mean','Median','Cumulative','0.5'},'Location','best','Box','off');
    end
end

function [x, y] = quest_to_xy(quest_results)
% 将 QUEST 结果转为拟合所需的样本点
% x 刺激强度
% y 正确性 1 表示正确 0 表示错误
    x=[]; y=[];
    for k=1:numel(quest_results)
        if isempty(quest_results(k).Stim_Actual), continue; end
        x(end+1,1)=quest_results(k).Stim_Actual;
        y(end+1,1)=double(quest_results(k).Correct);
    end
end

function fit = fit_psychometric(x, y, model)
% 心理物理函数拟合
% 使用最大似然目标函数对参数 alpha beta 求解
% gamma 固定为 0.5 表示二选一任务的猜测率
% lambda 固定为 0 表示不考虑失误率
    gamma=0.5; lambda=0;
    x=x(:); y=y(:);
    m=~isnan(x)&~isnan(y); x=x(m); y=y(m);

    if isempty(x)
        fit=struct('alpha',NaN,'beta',NaN,'nll',NaN,'model',model); return;
    end

    a0=max(median(x), eps);
    b0=max((max(x)-min(x))/5, 0.2);
    theta0=[log(a0), log(b0)];
    opt=optimset('Display','off');
    nll_fun=@(th) negll_psy(th,x,y,model,gamma,lambda);
    th_hat=fminsearch(nll_fun,theta0,opt);

    fit.alpha=exp(th_hat(1));
    fit.beta =exp(th_hat(2));
    fit.nll  =nll_fun(th_hat);
    fit.model=model;
end

function nll = negll_psy(th, x, y, model, gamma, lambda)
% 负对数似然
% p 为给定参数下的反应概率
% 为数值稳定将 p 限制在 1e-6 到 1-1e-6
    alpha=exp(th(1)); beta=exp(th(2));
    p=psychometric_p(x,alpha,beta,model);
    p=gamma+(1-gamma-lambda).*p;
    p=min(max(p,1e-6),1-1e-6);
    nll=-sum(y.*log(p)+(1-y).*log(1-p));
end

function p = psychometric_p(x, alpha, beta, model)
% 不同心理物理函数形式
% logistic 逻辑函数
% weibull  Weibull 累积分布形式
% normal   正态分布累积分布函数形式
    switch lower(model)
        case 'logistic', p = 1 ./ (1 + exp(-(x-alpha)./beta));
        case 'weibull',  p = 1 - exp(-(x./alpha).^beta);
        case 'normal',   p = 0.5*(1 + erf(((x-alpha)./beta)./sqrt(2)));
        otherwise, error('未知模型: %s', model);
    end
end

function p = Transfer(x, threshold)
% 任务正确率随刺激强度变化的转移函数
% threshold 越大表示阈值越高
% 输出范围 0.5 到 1
    p = 0.5 + 0.5*(1 - exp(-x/(2*threshold)));
end
