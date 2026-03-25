# KF_analysis

## 项目目标
项目的核心目标是把matlab的分析转换成python的版本，最后在vscode的jupyter文件中完成完整的版本。我希望新版本和matlab版本的代码实现的功能内核完全一致 但是更清晰更简洁。

test.ipynb是需要修改成的python版本

这个仓库里的 MATLAB 分析，核心上是“运动轨迹（左右手）→ X/Y 维度对齐与误差建模 → 群体统计/可视化”，并且在 quest 条件下单独走一条心理物理阈值估计路径。最合理的运行顺序是：先 main2.m（构建 results 主数据结构）、再 main_test.m（基于 results 做统计）、最后可选 delay_analysis.m（基于 data_import 做更细的延迟/误差方向分析），而 negLogLikelihoodr.m 是前两条路径都调用的底层似然函数。这个依赖关系分别体现在：main2.m 自己创建 results 并在末尾画图/单被试分析；main_test.m 直接读取 results；delay_analysis.m 开头直接写 data = data_import，说明它通常依赖 main2.m 先把数据加载到工作区。 

## main2.m（主流程，最重要）
main2.m 的目标非常完整：读每个文件、识别被试和条件、对轨迹做预处理、计算速度互相关（CCG）与时间滞后、用最大似然估计 Kalman 观测噪声标准差作为不确定性，再对 QUEST 做阈值后验与心理测量曲线拟合，最后输出单被试和群体可视化。也就是说，它不仅是“算指标”，还是“出图脚本”。

在数据层面，它先用 exam_load() 读取结构体数组，每个元素含 filename 和 c3d，再限定被试 101–116 和五个条件（vis.kinarm / low.kinarm / mid.kinarm / high.kinarm / quest.kinarm）。每个文件解析出 subject_id 与 condition 后，写到 results.subject_XXX.<short_cond> 层级字段里（vis/low/mid/hig/que）。

你特别关心的 X/Y 变量处理 在 main2.m 里是“并行但独立”的：

先把右手当 response（rx, ry），左手当 target（lx, ly），从试次 3 到 300 连起来，形成连续时序。

再组织为 2 × T 矩阵：第 1 行 X，第 2 行 Y（response_ori=[rx';ry'], target_ori=[lx';ly']）。

然后分别对 X 和 Y 计算速度互相关：先 diff 得速度、zscore 标准化，再 xcorr(...,'coeff') 得 ccg_x/ccg_y，各自取峰值索引得到 maxLag_x/maxLag_y。

再按 X 和 Y 自己的 lag 分别做时移对齐（lag_align），不是共用同一个 lag。

对齐后每 50 点抽 1 点（Fs=50），得到低频序列，再分别估计 r_std_x 与 r_std_y。

estimate_r_std 里，X 或 Y 单维序列都做同样处理：去均值、去掉开头 clip 段、换算到 cm 后，用 fminunc(@negLogLikelihoodr, ...) 优化 r_log，最终输出 sqrt(exp(r_log)) 作为该维度不确定性（cm）。这一步是你后续所有“不确定性统计”的源头。 

quest 条件则完全是另一条路径：按轮次读关键试次（第4个试次给刺激强度，第6个试次给选择），依据正确/错误对阈值后验做贝叶斯更新，最后从后验分布提取 MAP/mean/median，并把 trial 级 (stim, correct) 转成拟合数据，分别拟合 logistic / weibull / normal 三种心理物理函数的 alpha,beta 参数。你在 main_test.m 里看到的 6 类阈值来源就在这里。 

可视化方面，main2.m 一次性提供：群体 CCG（X/Y均值+SEM 阴影）、群体不确定性（X/Y散点+均值条+误差线）、群体 QUEST 指标散点，以及交互式单被试轨迹（1D X/Y 和 2D 轨迹渐变）+单被试汇总图（四条件 CCG 峰值、四条件 X/Y sigma、QUEST 后验图与拟合参数文本）。这说明它是你分析工程里的“主入口+主报告脚本”。 

## main_test.m（统计整合脚本，基于 results）
main_test.m 不再重新跑轨迹处理，而是把 main2.m 的 results 当输入，专门做统计学比较与汇总展示。结构上分两大块：第一块是“QUEST阈值 vs mid 条件 X 不确定性”的相关；第二块是“X/Y 不确定性、CCG峰值、lag”的组统计和显著性绘图。 

第一块里最关键是横轴定义：
x_plot = norminv(0.75,0,1) * mid_sigma_x，即把 mid.r_std_x 乘 75%分位常数换到 JND 尺度；纵轴则是六种阈值（后验MAP、后验均值、中位数、以及三种拟合 alpha）。每个阈值都和这条 X 轴做相关，画散点、线性拟合线、并标注被试编号。这个设计本质是在问：“中距离本体 X 方向噪声是否能预测 QUEST 阈值估计？” 

第二块中，collect_group_stats 统一抽取被试×条件矩阵：

不确定性：rstd_x, rstd_y

CCG峰值：ccgmax_x, ccgmax_y（直接取 ccg_x/ccg_y 最大值）

峰值 lag：lag_x, lag_y（由峰值索引 (idx-51)/50 转秒）
这就把所有 X/Y 指标整理到同一统计格式。 

随后的统计逻辑很系统：

近/中/远条件下，X 和 Y 分开做 RM-ANOVA；若主效应显著，再做近-中、近-远、中-远配对比较；

每个条件做 X vs Y 的配对 t（不确定性一套、CCG max 再一套）；

只比较 vis 与 mid 的 lag（X/Y 分开配对 t）；

统一转成分组柱图（均值±SE）并叠加显著性线/星号。

所以从功能定位看：main2.m 负责“算”，main_test.m 负责“检验与论文图式汇总”。

## delay_analysis.m（误差时间结构 + 方向分布拓展）
delay_analysis.m 的定位不是替代 main2.m，而是扩展“时间误差结构”和“误差方向偏置”分析。它依赖现成 data_import，重新遍历非 quest 条件并构建另一套 results/summary。 

X/Y 处理在这里更偏向“几何误差”而不是“互相关峰值统计”：

先拼接 rx,ry,lx,ly，再进入两个子模块：

do_kalman_segments：每 50 点下采样后按 20 点分段，X/Y 各自拟合 r_std_x/y，再取均值（并过滤 >10）；

do_ccg_align_and_distance：转 cm 后计算误差。

do_ccg_align_and_distance 是本文件最关键：

先设 tx,ty,rx,ry，并有一个额外 rx = rx - 50 的 X 偏移校正；

计算了 xcorr 但把 lag 强制设为 0（即最终不做时移补偿）；

定义有符号误差 ex = arx-atx, ey = ary-aty；

再派生绝对距离 dx=|ex|, dy=|ey|、欧式距离 d=sqrt(ex^2+ey^2)；

方向角 theta = atan2(ey, ex)，范围 (-pi,pi]；

输出三层时间表征：全程逐点（300s相关）、20s段均值、20s段内逐点平均。

然后它在 summary 聚合到群体层面，做每条件均值/标准差/SEM，最后出三类图：

单被试：300s全程误差（欧式/X/Y）+ 20s 段内平均；

群体：X/Y 20s 均值±SEM；

误差积累分析：取 15–16s 与 18–19s 的均值差，低/中/远做 ANOVA（X 与 Y 分开）；

方向分布：按条件画 error_theta 极坐标直方图，并用合成向量长度 R 表征方向一致性。

negLogLikelihoodr.m（底层函数，供 Kalman 拟合）
这个函数就是“给定 Q 和观测噪声参数 rr(log R)，计算目标序列 X 与估计序列 Xhat 在 Kalman 模型下的负对数似然”。内部先 rr=exp(rr)，算后验方差 pp、Kalman 增益 k 和差分矩阵 d，再累计每列数据的似然。main2.m 和 delay_analysis.m 都通过 fminunc 反复调用它来反推最优 R
