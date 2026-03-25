%% ============ 全局绘图风格（与main2一致：集中管理）============
STYLE = viz_style("science");
set_global_style(STYLE);

%% ============ 基本设置 ============
excluded_subjects = [110]; % 需要排除的被试ID（按需修改）
fprintf('\n正在进行相关分析，将排除以下被试: %s\n', join_num(excluded_subjects));

%% ============ Part 1：6种阈值 vs mid条件X方向不确定性相关 ============
[subjects, mid_sigma_x, thr, thr_name] = collect_quest_mid_x(results, excluded_subjects);
fprintf('相关分析包含的有效被试: %s\n', join_num(subjects));

x_plot = norminv(0.75, 0, 1) .* mid_sigma_x; % 横轴：σ -> JND尺度（保持原逻辑）
fig1 = figure('Position', [100, 100, 1500, 820], 'Color', 'w', 'Name', 'QUEST阈值相关');
tl1 = tiledlayout(fig1, 2, 3, 'TileSpacing','compact', 'Padding','compact');

fprintf('\n========== 6种阈值与mid \\sigma_x 的相关结果 ==========\n');
for k = 1:numel(thr)
    ax = nexttile(tl1, k);
    yk = thr{k};

    ok = isfinite(x_plot) & isfinite(yk);
    x = x_plot(ok); y = yk(ok); sid = subjects(ok);

    scatter(ax, x, y, 46, 'filled', ...
        'MarkerFaceColor', STYLE.C(1,:), 'MarkerFaceAlpha', 0.85, ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
    hold(ax, 'on');

    if numel(x) >= 3
        [r, p] = corr(x, y);
        r2 = r^2;

        pp = polyfit(x, y, 1);
        xx = linspace(min(x)*0.95, max(x)*1.05, 100);
        plot(ax, xx, polyval(pp, xx), '-', 'Color', STYLE.C(4,:), 'LineWidth', 1.4);

        % 被试编号：小字号+深灰，避免喧宾夺主
        for i = 1:numel(x)
            text(ax, x(i), y(i), num2str(sid(i)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                'FontSize', 9, 'Color', [0.25 0.25 0.25]);
        end

        title(ax, sprintf('%s\nr=%.3f, p=%.3g, R^2=%.3f', thr_name{k}, r, p, r2));
        fprintf('%s: r=%.3f, p=%.3g, R^2=%.3f, n=%d\n', thr_name{k}, r, p, r2, numel(x));
    else
        title(ax, sprintf('%s（有效点不足）', thr_name{k}));
        fprintf('%s: 有效点不足(n=%d)\n', thr_name{k}, numel(x));
    end

    xlabel(ax, 'mid 条件下 X 方向不确定性（JND尺度）');
    ylabel(ax, '阈值 (cm)');
    apply_axes_style(ax, STYLE);
end
sgtitle(tl1, '6种QUEST阈值与mid条件下X方向不确定性的相关分析', 'FontWeight','bold');

%% ============ Part 2：统计分析（不确定性、CCG、lag） ============
fig2 = figure('Position', [100, 100, 1250, 760], 'Color','w', 'Name', '统计分析结果');
tl2 = tiledlayout(fig2, 2, 2, 'TileSpacing','compact', 'Padding','compact');

cond_names  = {'vis', 'low', 'mid', 'hig'};
cond_labels = {'视觉', '近',  '中',  '远'};

[sub_all, rstd_x, rstd_y, ccgmax_x, ccgmax_y, lag_x, lag_y] = collect_group_stats(results, cond_names);

%% ------- 分析1：近/中/远下不确定性差异（X/Y分别RM-ANOVA）-------
idx_complete = all(isfinite(rstd_x(:,2:4)),2) & all(isfinite(rstd_y(:,2:4)),2);
complete_data_available = nnz(idx_complete) >= 3;

sig_x = []; sig_y = [];
p_x = NaN; p_y = NaN; tbl_x = {}; tbl_y = {};

if complete_data_available
    x_lmh = rstd_x(idx_complete, 2:4);
    y_lmh = rstd_y(idx_complete, 2:4);

    [p_x, tbl_x] = anova_rm(x_lmh, 'off'); % X方向RM-ANOVA
    [p_y, tbl_y] = anova_rm(y_lmh, 'off'); % Y方向RM-ANOVA

    if p_x(1) < 0.05, sig_x = posthoc_pairs(x_lmh); end
    if p_y(1) < 0.05, sig_y = posthoc_pairs(y_lmh); end
else
    warning('没有足够的被试有完整的近/中/远数据，无法进行RM-ANOVA。');
end

%% ------- 分析2：各条件下 X vs Y 不确定性（配对t）-------
paired_xy = paired_t_allconds(rstd_x, rstd_y, cond_names, cond_labels);
fprintf('\n\n========== 分析2：不同条件下X vs Y方向不确定性比较 ==========\n');
print_paired_xy(paired_xy);

%% ------- 分析4：vis vs mid 的 lag（X/Y分别配对t）-------
[lag_stat_x, lag_stat_y] = compare_vis_mid(lag_x, lag_y);

%% ------- 分析5+6：各条件下 CCG max 的 X vs Y（配对t）-------
ccg_xy = paired_t_allconds(ccgmax_x, ccgmax_y, cond_names, cond_labels);
fprintf('\n\n========== 分析5+6：不同条件下CCG最高相关值比较（X vs Y） ==========\n');
print_paired_xy(ccg_xy);

%% ============ 绘图A：各条件 X/Y 不确定性（均值±SE + 显著性） ============
axA = nexttile(tl2, [1 2]);
order   = [2 3 4 1];                 % 近/中/远/视觉
labelsA = {'近','中','远','视觉'};

[meanX, seX, pXY] = pack_bar(paired_xy, cond_names, order, 'x');
[meanY, seY, ~  ] = pack_bar(paired_xy, cond_names, order, 'y');

bar_group_se(axA, meanX, seX, meanY, seY, labelsA, STYLE, '不确定性 (σ)');
title(axA, '不同条件下X与Y方向本体感觉不确定性比较', 'FontWeight','bold');

add_sig_xy(axA, meanX, seX, meanY, seY, pXY); % 每个条件的X-Y显著性

% 近/中/远：RM-ANOVA后事后比较（只画显著）
if complete_data_available && ~isempty(sig_x)
    add_sig_pairs(axA, sig_x, -0.18, STYLE.C(4,:)); % X组(左)
end
if complete_data_available && ~isempty(sig_y)
    add_sig_pairs(axA, sig_y, +0.18, STYLE.C(1,:)); % Y组(右)
end

%% ============ 绘图B：vis vs mid lag（均值±SE + 显著性） ============
axB = nexttile(tl2, 3);
bar_group_se(axB, lag_stat_x.mean, lag_stat_x.se, lag_stat_y.mean, lag_stat_y.se, ...
    {'视觉','本体(中)'}, STYLE, '时间延迟 (s)');
title(axB, '视觉 vs 本体(中) 条件下CCG最高点时间比较', 'FontWeight','bold');
add_sig_two(axB, lag_stat_x.p, -0.18, lag_stat_x.mean, lag_stat_x.se);
add_sig_two(axB, lag_stat_y.p, +0.18, lag_stat_y.mean, lag_stat_y.se);

%% ============ 绘图C：各条件 CCG max（均值±SE + 显著性） ============
axC = nexttile(tl2, 4);
[ccgMeanX, ccgSeX, ccgPXY] = pack_bar(ccg_xy, cond_names, 1:4, 'x');
[ccgMeanY, ccgSeY, ~     ] = pack_bar(ccg_xy, cond_names, 1:4, 'y');

bar_group_se(axC, ccgMeanX, ccgSeX, ccgMeanY, ccgSeY, {'视觉','近','中','远'}, STYLE, '相关系数');
title(axC, '不同条件下X vs Y方向CCG最高相关值比较', 'FontWeight','bold');
ylim(axC, [0 1.1]);
add_sig_xy(axC, ccgMeanX, ccgSeX, ccgMeanY, ccgSeY, ccgPXY);

sgtitle(tl2, '本体感觉不确定性与跨相关分析', 'FontWeight','bold');

%% ============ 打印ANOVA摘要 ============
fprintf('\n\n========== 分析1：近/中/远下卡尔曼不确定性RM-ANOVA ==========\n');
if ~complete_data_available
    fprintf('数据不足：无法进行RM-ANOVA。\n');
else
    fprintf('X方向: F(%d,%d)=%.2f, p=%.4f\n', tbl_x{3,3}, tbl_x{4,3}, tbl_x{3,5}, p_x(1));
    if ~isempty(sig_x), print_sig_pairs(sig_x, {'近','中','远'}); end

    fprintf('Y方向: F(%d,%d)=%.2f, p=%.4f\n', tbl_y{3,3}, tbl_y{4,3}, tbl_y{3,5}, p_y(1));
    if ~isempty(sig_y), print_sig_pairs(sig_y, {'近','中','远'}); end
end
fprintf('\n********* 统计分析完成 *********\n');

%% ======================================================================
%                               函数区
%% ======================================================================

function [subjects, mid_sigma_x, thr, thr_name] = collect_quest_mid_x(results, excluded_subjects)
% 提取：mid条件X方向不确定性 + 6种QUEST阈值
subjects = []; mid_sigma_x = [];
thr_map = []; thr_mean = []; thr_median = [];
thr_fit_logistic = []; thr_fit_weibull = []; thr_fit_normal = [];

fields = fieldnames(results);
for f = 1:numel(fields)
    sf = fields{f};
    if ~startsWith(sf, 'subject_'), continue; end
    sid = str2double(erase(sf, 'subject_'));
    if ismember(sid, excluded_subjects)
        fprintf('被试 %d 已从相关分析中排除\n', sid);
        continue;
    end

    if ~isfield(results.(sf), 'mid') || ~isfield(results.(sf), 'que'), continue; end
    if ~isfield(results.(sf).mid, 'r_std_x'), continue; end

    q = results.(sf).que;
    need = {'quest_results','threshold_val','threshold_mean','threshold_median','fit_logistic','fit_weibull','fit_normal'};
    if ~all(isfield(q, need)), continue; end

    final_prob = q.quest_results(end).Threshold_Prob;
    final_prob = final_prob ./ sum(final_prob);
    [~, map_idx] = max(final_prob);
    threshold_map = q.threshold_val(map_idx);

    subjects(end+1,1) = sid;
    mid_sigma_x(end+1,1) = results.(sf).mid.r_std_x;

    thr_map(end+1,1) = threshold_map;
    thr_mean(end+1,1) = q.threshold_mean;
    thr_median(end+1,1) = q.threshold_median;

    thr_fit_logistic(end+1,1) = q.fit_logistic.alpha;
    thr_fit_weibull(end+1,1)  = q.fit_weibull.alpha;
    thr_fit_normal(end+1,1)   = q.fit_normal.alpha;
end

thr = {thr_map, thr_mean, thr_median, thr_fit_logistic, thr_fit_weibull, thr_fit_normal};
thr_name = {'后验MAP', '后验均值', '后验中位数', 'Logistic拟合', 'Weibull拟合', 'Normal拟合'};
end

function [subjects_all, rstd_x, rstd_y, ccgmax_x, ccgmax_y, lag_x, lag_y] = collect_group_stats(results, cond_names)
% 提取：不确定性(x/y)、CCG max(x/y)、lag(x/y)
sf = fieldnames(results);
sf = sf(startsWith(sf,'subject_'));
subjects_all = unique(str2double(erase(sf,'subject_')));

nS = numel(subjects_all); nC = numel(cond_names);
rstd_x = nan(nS,nC); rstd_y = nan(nS,nC);
ccgmax_x = nan(nS,nC); ccgmax_y = nan(nS,nC);
lag_x = nan(nS,nC); lag_y = nan(nS,nC);

for s = 1:nS
    key = ['subject_' num2str(subjects_all(s))];
    if ~isfield(results, key), continue; end

    for c = 1:nC
        cn = cond_names{c};
        if ~isfield(results.(key), cn), continue; end

        rstd_x(s,c) = getfield_nan(results.(key).(cn), 'r_std_x');
        rstd_y(s,c) = getfield_nan(results.(key).(cn), 'r_std_y');

        ccgx = getfield_vec(results.(key).(cn), 'ccg_x');
        ccgy = getfield_vec(results.(key).(cn), 'ccg_y');

        if ~isempty(ccgx) && all(isfinite(ccgx))
            [ccgmax_x(s,c), idx] = max(ccgx);
            lag_x(s,c) = (idx - 51) / 50; % 50Hz，中心点51（保持原逻辑）
        end
        if ~isempty(ccgy) && all(isfinite(ccgy))
            [ccgmax_y(s,c), idx] = max(ccgy);
            lag_y(s,c) = (idx - 51) / 50;
        end
    end
end
end

function v = getfield_nan(s, name)
if isstruct(s) && isfield(s, name) && isfinite(s.(name)), v = s.(name); else, v = NaN; end
end

function v = getfield_vec(s, name)
if isstruct(s) && isfield(s, name), v = s.(name); else, v = []; end
end

function out = join_num(x)
if isempty(x), out = '(none)'; return; end
out = strjoin(arrayfun(@num2str, x(:)', 'UniformOutput', false), ', ');
end

function [p, table, stats] = anova_rm(y, displayopt)
% 重复测量ANOVA：行=被试，列=条件
if nargin < 2, displayopt = 'on'; end
[n, k] = size(y);

y_mean = mean(y);
subj_mean = mean(y, 2);
grand_mean = mean(y(:));

SS_total = sum((y(:) - grand_mean).^2);
SS_subj  = sum((subj_mean - grand_mean).^2) * k;
SS_treat = sum((y_mean - grand_mean).^2) * n;
SS_error = SS_total - SS_subj - SS_treat;

df_subj  = n - 1;
df_treat = k - 1;
df_error = (n - 1) * (k - 1);

MS_treat = SS_treat / df_treat;
MS_error = SS_error / df_error;

F_treat = MS_treat / MS_error;
p_treat = 1 - fcdf(F_treat, df_treat, df_error);

table = cell(4, 6);
table(1,:) = {'Source','SS','df','MS','F','p'};
table(2,:) = {'Subjects',  SS_subj,  df_subj,  SS_subj/df_subj,  '',     ''};
table(3,:) = {'Treatments',SS_treat, df_treat, MS_treat,       F_treat, p_treat};
table(4,:) = {'Error',     SS_error, df_error, MS_error,       '',     ''};

p = [p_treat; NaN; NaN];
stats = struct('n', n, 'k', k, 'means', y_mean(:)', 'df', df_error, 's', sqrt(MS_error));
if strcmp(displayopt, 'on'), disp('重复测量ANOVA表:'); disp(table); end
end

function sig = posthoc_pairs(M)
% 近-中、近-远、中-远：配对t检验，返回 [i j p]
[~, p12] = ttest(M(:,1), M(:,2));
[~, p13] = ttest(M(:,1), M(:,3));
[~, p23] = ttest(M(:,2), M(:,3));
sig = [1 2 p12; 1 3 p13; 2 3 p23];
end

function res = paired_t_allconds(Ax, Ay, cond_names, cond_labels)
% 每个条件：Ax vs Ay 配对t检验
res = struct();
for c = 1:numel(cond_names)
    cn = cond_names{c}; lab = cond_labels{c};
    ok = isfinite(Ax(:,c)) & isfinite(Ay(:,c));
    if nnz(ok) < 3
        fprintf('\n%s条件下数据不足，无法进行配对t检验\n', lab);
        res.(cn).valid = false; res.(cn).label = lab;
        continue;
    end

    x = Ax(ok,c); y = Ay(ok,c);
    [~, p, ci, stats] = ttest(x, y);
    d = mean(x - y) / std(x - y);

    res.(cn) = struct('valid',true,'label',lab,'t',stats.tstat,'df',stats.df,'p',p,'ci',ci,'d',d, ...
        'x_mean',mean(x),'y_mean',mean(y),'x_std',std(x),'y_std',std(y),'n',nnz(ok));
end
end

function print_paired_xy(res)
conds = fieldnames(res);
for i = 1:numel(conds)
    cn = conds{i};
    if ~res.(cn).valid
        fprintf('\n%s: 数据不足\n', res.(cn).label);
        continue;
    end
    R = res.(cn);
    fprintf('\n%s条件 (X vs Y): t(%d)=%.2f, p=%.4f, d=%.2f, n=%d\n', R.label, R.df, R.t, R.p, R.d, R.n);
    fprintf('  均值差(X-Y)=%.3f, 95%%CI=[%.3f, %.3f]\n', R.x_mean-R.y_mean, R.ci(1), R.ci(2));
    fprintf('  X均值(SD)=%.3f(%.3f), Y均值(SD)=%.3f(%.3f)\n', R.x_mean, R.x_std, R.y_mean, R.y_std);
end
end

function [meanV, seV, pV] = pack_bar(res, cond_names, order, which)
% 按顺序打包bar数据（均值、SE、p）
meanV = nan(1,numel(order)); seV = meanV; pV = meanV;
for i = 1:numel(order)
    cn = cond_names{order(i)};
    if ~isfield(res, cn) || ~res.(cn).valid, continue; end
    R = res.(cn);
    if which == 'x'
        meanV(i) = R.x_mean; seV(i) = R.x_std / sqrt(R.n);
    else
        meanV(i) = R.y_mean; seV(i) = R.y_std / sqrt(R.n);
    end
    pV(i) = R.p;
end
end

function bar_group_se(ax, meanX, seX, meanY, seY, xticklabels, STYLE, ylab)
% 分组bar + 误差线：稳定兼容（不依赖 ax.Children 顺序）
data = [meanX(:), meanY(:)];
h = bar(ax, data, 'grouped'); hold(ax,'on');

% 颜色与透明度：X=橙，Y=蓝（与其它文件一致）
if numel(h) >= 2
    h(1).FaceColor = STYLE.C(4,:); h(1).EdgeColor = 'none'; h(1).FaceAlpha = 0.85;
    h(2).FaceColor = STYLE.C(1,:); h(2).EdgeColor = 'none'; h(2).FaceAlpha = 0.85;
end

% 误差线位置：优先用 XEndPoints（新版本），否则手动估算（旧版本）
[x1, x2] = bar_centers(h, numel(meanX));
errorbar(ax, x1, meanX, seX, 'k', 'LineStyle','none', 'LineWidth',1.1);
errorbar(ax, x2, meanY, seY, 'k', 'LineStyle','none', 'LineWidth',1.1);

set(ax, 'XTick', 1:numel(meanX), 'XTickLabel', xticklabels);
ylabel(ax, ylab);
legend(ax, {'X方向','Y方向'}, 'Location','best');
apply_axes_style(ax, STYLE);

m = max([meanX(:)+seX(:); meanY(:)+seY(:)], [], 'omitnan');
ylim(ax, [0, m*1.35 + eps]); % 给显著性标记留空间
end

function [x1, x2] = bar_centers(h, n)
% 返回两组bar在每个类别上的中心x坐标（兼容新旧版本）
try
    x1 = h(1).XEndPoints(:);
    x2 = h(2).XEndPoints(:);
catch
    % 旧版：手动估算 grouped bar 的偏移
    groupwidth = min(0.8, 2/(2+1.5));
    x = (1:n)';
    x1 = x - groupwidth/4;
    x2 = x + groupwidth/4;
end
end

function add_sig_xy(ax, meanX, seX, meanY, seY, pvals)
% 每个条件：画 X vs Y 的显著性（线 + 星号/ns）
yr = range(ylim(ax));
for i = 1:numel(pvals)
    if ~isfinite(pvals(i)), continue; end
    y0 = max(meanX(i)+seX(i), meanY(i)+seY(i)) + 0.05*yr;
    plot(ax, [i-0.18 i+0.18], [y0 y0], 'k', 'LineWidth', 1);
    text(ax, i, y0 + 0.01*yr, pstar(pvals(i)), 'HorizontalAlignment','center', 'FontSize', 12);
end
end

function add_sig_pairs(ax, sig_mat, xshift, col)
% 近/中/远之间：只画显著的事后比较（同一组内部）
yr = range(ylim(ax));
ybase = max(ylim(ax)) - 0.18*yr;
h = 0;
for s = 1:size(sig_mat,1)
    p = sig_mat(s,3);
    if ~isfinite(p) || p >= 0.05, continue; end
    i = sig_mat(s,1); j = sig_mat(s,2);
    y = ybase + h*0.06*yr;
    plot(ax, [i+xshift j+xshift], [y y], 'Color', col, 'LineWidth', 1.1);
    text(ax, (i+j)/2 + xshift, y + 0.01*yr, pstar(p), 'HorizontalAlignment','center', 'FontSize', 11, 'Color', col);
    h = h + 1;
end
end

function add_sig_two(ax, p, xshift, means, ses)
% 两条件：视觉 vs 本体中（在对应组上方画线+标记）
if ~isfinite(p), return; end
yr = range(ylim(ax));
y0 = max(means+ses) + 0.08*yr;
plot(ax, [1 2]+xshift, [y0 y0], 'k', 'LineWidth', 1.1);
text(ax, 1.5+xshift, y0 + 0.01*yr, pstar(p), 'HorizontalAlignment','center', 'FontSize', 12);
end

function s = pstar(p)
if p < 0.001, s = '***';
elseif p < 0.01, s = '**';
elseif p < 0.05, s = '*';
else, s = 'ns';
end
end

function print_sig_pairs(sig_mat, names)
fprintf('事后检验（配对t检验）:\n');
for i = 1:size(sig_mat,1)
    a = sig_mat(i,1); b = sig_mat(i,2); p = sig_mat(i,3);
    fprintf('  %s vs %s: p=%.4f%s\n', names{a}, names{b}, p, ternary(p<0.05,' (显著)',' (不显著)'));
end
end

function [sx, sy] = compare_vis_mid(lag_x, lag_y)
% vis(1) vs mid(3)：lag差异（X/Y分别配对t）
sx = lag_pair(lag_x(:,1), lag_x(:,3));
sy = lag_pair(lag_y(:,1), lag_y(:,3));

fprintf('\n\n========== 分析4：视觉 vs 本体(中) 条件下CCG最高点时间差异 ==========\n');
if sx.valid
    fprintf('X方向: vis=%.4f(SE=%.4f), mid=%.4f(SE=%.4f), p=%.4f, d=%.2f, n=%d\n', ...
        sx.mean(1), sx.se(1), sx.mean(2), sx.se(2), sx.p, sx.d, sx.n);
else
    fprintf('X方向：数据不足\n');
end
if sy.valid
    fprintf('Y方向: vis=%.4f(SE=%.4f), mid=%.4f(SE=%.4f), p=%.4f, d=%.2f, n=%d\n', ...
        sy.mean(1), sy.se(1), sy.mean(2), sy.se(2), sy.p, sy.d, sy.n);
else
    fprintf('Y方向：数据不足\n');
end
end

function S = lag_pair(vis_v, mid_v)
ok = isfinite(vis_v) & isfinite(mid_v);
if nnz(ok) < 3
    S = struct('valid',false,'mean',[NaN NaN],'se',[NaN NaN],'p',NaN,'d',NaN,'n',nnz(ok));
    return;
end
a = vis_v(ok); b = mid_v(ok);
[~, p] = ttest(a, b);
d = mean(a-b) / std(a-b);
S = struct('valid',true, 'mean',[mean(a) mean(b)], ...
    'se',[std(a)/sqrt(nnz(ok)) std(b)/sqrt(nnz(ok))], 'p',p, 'd',d, 'n',nnz(ok));
end

function t = ternary(cond, a, b)
if cond, t = a; else, t = b; end
end

function STYLE = viz_style(preset)
% 低饱和科研配色 + 全局参数（跨文件统一）
if nargin<1, preset="science"; end
STYLE.font = "Arial";
STYLE.font_size = 11;
STYLE.lw = 1.6;
STYLE.ax_lw = 1.0;
STYLE.grid_alpha = 0.12;

switch lower(string(preset))
    case "science"
        STYLE.C = [0.20 0.40 0.70;   % 蓝
                   0.85 0.37 0.01;   % 橙
                   0.18 0.63 0.17;   % 绿
                   0.55 0.35 0.64;   % 紫
                   0.35 0.35 0.35];  % 灰
    otherwise
        STYLE.C = [0.12 0.47 0.71; 0.20 0.63 0.17; 0.80 0.33 0.10; 0.60 0.31 0.64; 0.35 0.35 0.35];
end
end

function set_global_style(STYLE)
% 全局默认绘图样式：保证多个figure/子图一致
set(groot,"defaultFigureColor","w");
set(groot,"defaultAxesFontName",STYLE.font,"defaultAxesFontSize",STYLE.font_size);
set(groot,"defaultLineLineWidth",STYLE.lw);
set(groot,"defaultAxesLineWidth",STYLE.ax_lw);
set(groot,"defaultAxesTickDir","out");
set(groot,"defaultAxesBox","off");
end

function apply_axes_style(ax, STYLE)
% 坐标轴统一：淡网格、外刻度、去上右框
grid(ax, "on"); ax.GridAlpha = STYLE.grid_alpha;
ax.TickDir = "out"; ax.Box = "off";
end
