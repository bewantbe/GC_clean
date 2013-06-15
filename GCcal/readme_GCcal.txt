Granger因果关系分析 matlab 程序(精简版)说明文档

若要在别的目录使用这些代码(函数)只需在 Matlab 先执行 addpath <到目录 GCcal 的路径>, 这一步
也可在当前目录的 startup.m (对于 Octave 是 .octaverc) 中自动完成.
prj_neuron_gc/analyse_neuron.m 是示例文件(内含大量注释).
神经元数据生成程序的说明文档另见 readme_raster_tuning.txt

-------------------------------------------------------------------------------
数据结构:

p 变量数; od (或m) 阶数; len 数据长度; fftlen 傅里叶变换长度.

被分析序列 X 的格式:
        设数据时间长度是 len, 则 X 为 p 行 len 列矩阵(每列一个时刻).
系数矩阵列 A 的格式:
        A 为 p 行, p*od 列矩阵. (每个 p*p 的块是一个系数.
        从左起依次是 t-1, t-2, ... t-od 时刻的系数).
协方差序列 R 的格式:
        p行 p*(od+1) 列矩阵, od 是被分析的阶数
        左数第 1 块p*p方阵表示方差矩阵, 第 2 块表示时间差为 1 的协方差矩阵.
因果关系矩阵的格式:
        第 j 列第 i 行表示 j 变量对 i 变量的影响强度. 对角线恒是 0.
谱 S 的格式:
        p*p*fftlen 的 3-dim 数组

-------------------------------------------------------------------------------
函数功能 (方括号[]中的参数可省略):

gdata.m:          X = gdata(len, od, p[, noisecov, A])
        给定 序列长度, 阶数, 变量数, 参数矩阵, 噪音项方差矩阵 依据 AR 模型生成数据
        若省略后两参数, 则按内设的模型生成数据. 生成的开头 200 个数据会被丢弃.

para2cov.m        R = para2cov(od, p[, extm])
        给定 阶数, 变量数[, 额外阶数], 计算协方差序列 R. (使用了 fft)

getcov.m          R = getcov(X, od)
        依据数据 X 计算协方差序列. od 是需要计算的最大时间偏移量.
        即最多算到 cov(X(t),X(t-od))

getcovpd.m        R = getcovpd(X, od)
getcovpd2.m       R = getcovpd2(X, od)
        依据数据 X 计算协方差序列, 尽可能使总体协方差正定

nGrangerT.m       [Gc, Deps, Aall] = nGrangerT(X, od)
        依据给定的数据 X 及 想要分析的阶数 od 计算因果关系矩阵. 计算时减去了样本均值.
        返回值 Gc 是因果关系矩阵, Deps 是整体回归得到的方差矩阵. Aall 是拟合出的系数.

pos_nGrangerT.m           [Gc, Deps, Aall] = pos_nGrangerT(X, od)
pos_nGrangerT2.m   [Gc, Deps, Aall] = pos_nGrangerT2(X, od)
        功能同 nGrangerT, 但在数学上是完全正定的版本. 用于协助验证 nGrangerT 的计算结果的有效性.
        其中 pos_nGrangerT2 的时间和内存消耗与 nGrangerT 相仿(远小于pos_nGrangerT).
        若要显示条件数, 使用 pos_nGrangerT(X, od, 1) 或 pos_nGrangerT2(X, od, 1)

RGrangerT.m        [Gc, Deps, Aall] = RGrangerT(R)
        依据给定的协方差序列 R, 计算因果关系矩阵, 阶数由 R 的长短而定. 通常用于计算理论的 Gc 值,
        或是用于已知协方差求不同阶数回归的 Gc 值(比直接用 nGrangerT 节省计算量).

pairGrangerT.m     GC = pairGrangerT(X, m)
pairRGrangerT.m    GC = pairRGrangerT(R)
        两两计算 GC (非条件). pairGrangerT 用法同 nGrangerT, pairRGrangerT 用法同 RGrangerT.
        可用于比较条件 GC 与 非条件 GC 计算结果的差异.

gc_prob_nonzero.m   p = gc_prob_nonzero(gc, od, len)
gc_prob_intv.m      [gc_lower, gc_upper] = gc_prob_intv(gc, od, len[, a])
        计算 gc 值非零的概率(gc_prob_nonzero)和置信区间(gc_prob_intv).
        与 nGrangerT() 系列函数配套使用. 即是说假设自回归与联合回归所用的阶数相同.
        a 是置信水平, 默认是 0.95. 对于置信区间, 使用的是近似算法, 可能有百分之几的误差.
        gc_prob_nonzero 用 0~1 的数表示非零的概率, 1 表示 1 概率非零.
        输入的 gc 可以是 GC 矩阵也可以是一个 GC 分量.

chooseOrder.m      [best_od, xic] = chooseOrder(X[, ic_mode[, od_max]])
chooseROrder.m     [best_od, xic] = chooseROrder(R, len[, ic_mode])
chooseOrderFull.m  [od_joint, od_vec] = chooseOrderFull(X[, ic_mode[, od_max]])
        依据某些准则给出最佳阶数. ic_mode 是准则模式, 可以是 'BIC'(chooseOrder默认),
        'AIC'(chooseOrderFull默认) 或 'AICc'.
        chooseGCOrder, chooseGCOrderFull最大阶数由od_max指定, 默认设在99阶.
        od_max 直接影响了计算时间. chooseRGCOrder最大阶数则由协方差 R 的项数决定.
        可能没法识别<2的阶数.
        chooseGCOrderFull 同时给出联合回归所需的阶数od_joint, 及自回归(缺一回归)所需的阶数.
        对于坏条件的问题, 本组程序可能会给出错误结果. 参见函数 AnalyseSeries()

A2S.m              S = A2S(A, noisecov, fftlen)
        由系数 A 以及误差方差 noisecov 计算谱, 使用的 Fourier 变换长度是 fftlen.
        谱 S 是 p*p*fftlen 的三维数组.

S2cov.m            R = S2cov(S, od)
        由谱计算协方差. 需要至少保证 size(S,3) >= 2*(od+1), 一般取 size(S,3) > 8*od
        与 A2S 配合使用可以实现 AR 回归系数到协方差 R 的转换.

ARroots.m          la = ARroots(A)
        由系数 A 给出 AR 过程的所有特征根 la, 用于判断平稳性.

ARregression.m     [Aall, Deps] = ARregression(R)
        计算由协方差 R 最小二乘回归得到的系数 Aall 和 残差 Deps. 阶数由 R 的长度确定.
        用于不需要计算 GC 的场合.

AnalyseSeries.m    [oGC, oDe, R] = AnalyseSeries(X, s_od[, bad_mode])
        计算用不同阶数拟合得到的 GC, 以及误差方差 De, 所需的协方差 R.
        oGC, oDe 都是 p*p*length(s_od) 阶三维数组.
        若问题是坏条件的, 可以尝试以 bad_mode = 1 调用. bad_mode 默认是 0. 
        若 bad_mode = 1,  R 仍是原来算法得到的结果.

AnalyseSeries2.m   [aic_od, bic_od, zero_GC, oAIC, oBIC] = AnalyseSeries2(s_od, oGC, oDe, len)
        给出进一步的阶数分析结果. 配合 AnalyseSeries.m 使用. aic_od 是AIC阶数, bic_od 是BIC阶数;
        zero_GC 是插值到阶数为零时的GC. oAIC, oBIC 是各个阶数的 AIC, BIC 值.

GC_regression.m    (Matlab script)
        对 X 执行 AnalyseSeries 和 AnalyseSeries2 并求其自回归残差 srd, 联合回归残差 rd.
        
WhiteningFilter.m  [srd, sas] = WhiteningFilter(X, use_od)
        返回 X 的自回归噪音项(srd), 也即"白化" X. sas是白化使用的系数. 自回归阶数统一为 use_od.

nGrangerF.m         wGC = nGrangerF(X[, od[, fftlen]])
        由原始数据计算序列两两间的频域条件 GC, 支持多变量. od 是内部回归时用的阶数,
        fftlen 是频域长度. 内部调用 singleRGrangerF.m 进行计算.

singleGrangerF.m    wGc = singleGrangerF(X, id_y, id_x[, fftlen])
singleRGrangerF.m   wGc = singleRGrangerF(R, id_y, id_x[, od[, fftlen]])
        计算群组 id_y (对应 X(id_y,:)) 到群组 id_x (对应 X(id_x,:)) 的频域GC.
        向量 id_x, id_y 存的是属于该群组的数据序列编号, 从 1 开始.
        singleRGrangerF 仅使用协方差 R 数据.
        od 表示使用的阶数, 可以分别指定联合回归与自回归的阶数(如od=[3 7]).

nGrangerF2.m
[gcy2x, gcx2y, fx2y, fy2x, Sxx, Syy] = nGrangerF2(X, od, fftlen)
        由原始数据计算频域 GC, 仅支持两变量. od 是内部回归时用的阶数. 回归后调用 nGrangerFA.m

nGrangerFA2.m
[gcy2x, gcx2y, fx2y, fy2x, Sxx, Syy] = nGrangerFA2(A, noisecov, fftlen)
        由联合回归的系数 A 以及误差方差 noisecov, 计算频域 GC, 只用于二维.
        返回值分别是 GC(y->x), GC(x->y), 两个频域 GC 的分量, 两个序列的自谱.

--------------------------------------
以下是用于实验的文件, 可靠性不保证, 慎用

experimental_tools/gen_hfreq_coef.m    [G, de] = gen_hfreq_coef(r, f[, od])
        生成具有(几乎)单一频率 f 的 AR 系数 G. de 是误差项, 使得结果的方差是 1.

experimental_tools/partial_GC.m        pGC = partial_GC(srdx, srdy, od)
        ?偏相关GC? 只用 y (而不用x)的历史预测 x_t 得到 y 对 x 的"因果关系".

experimental_tools/conv1mat.m          C = conv1mat(A, B)
        计算矩阵序列的卷积. A,B,C 都是平铺的矩阵列.

experimental_tools/X2Sxx.m             z = X2Sxx(X, fftlen)
        由数据直接计算自谱.

experimental_tools/shuffleData.m       X = shuffleData(X)
        打乱各数据点的时间.

experimental_tools/hilbert.m           f=hilbert(f, N = [], dim = [])
        返回实信号的复表示, f 的虚部是 f 的希尔伯特变换, 实部是原始信号本身.
        abs(hilbert(f)) 则得到信号的包络线.
        Matlab 中已有该函数

如果遇到分析出的因果关系出现复数或小于零或NaN, 可能问题是坏条件的, 请使用 pos_nGrangerT.m 或
pos_nGrangerT2.m 以确保解的是正定问题, 并且顺便查看矩阵的条件数. 另外对数据作差分也常常能显著
减小问题的条件数(需要的阶数似乎会升高).
但这时问题本质上仍是坏条件的, 所得结果仍需谨慎看待.


===============================================================================
使用示例:

% 添加 GCcal 为工作目录, 这一步也可在 startup.m (对于 Octave 是 .octaverc) 中自动完成, 下同
% addpath GCcal                       

% 计算长 10000 的序列, [Ding] 中"经典" 三元二阶模型的 Granger Causality
nGrangerT(gdata(10000,2,3), 2)

% 计算上面结果的理论值(序列长度趋于无穷)
RGrangerT(para2cov(2,3))

% 为了直观的比较, 显示多几位精度
format long

% 计算"无穷阶"的结果(即是谱的结果). 其中30是任取的一个较大阶数(实际上计算的是 30+2 阶)
RGrangerT(para2cov(2,3,30))

% 做上面的计算的另一种方式
noisecov = diag([0.3 1.0 0.2]);
A = [-0.8  0.0 -0.4  0.5 -0.2  0.0;
      0.0 -0.9  0.0  0.0  0.8  0.0;
      0.0 -0.5 -0.5  0.0  0.0  0.2];
S = A2S(A, noisecov, 1024);           % 算谱
RGrangerT(S2cov(S, 2))                % 由谱得到协方差再到 GC
RGrangerT(S2cov(S, 32))


===============================================================================
% 计算神经元的 GC
% 计算/分析神经元数据的完整例子见 prj_neuron_gc/analyse_neuron.m
% (设工作目录为 prj_neuron_gc/ )
addpath ../GCcal                               % 添加 GCcal 为工作目录
mt = [0 0 0; 1 0 0; 0 1 0];                    % 1->2->3
save('-ascii', 'cormat_tmp.txt', 'mt');        % 保存连接关系矩阵到文本文件 cormat_tmp.txt
% 开始计算. 不显示图形, 显示即时运行信息, 3神经, 时间10000ms, 读取矩阵文件 cormat_tmp.txt
system('./raster_tuning -ng -v -n 3 -t 10000 -mat cormat_tmp.txt');
X = load('data/staffsave.txt')';               % 注意转置
nGrangerT(X, 10)


===============================================================================
% 计算某两群组间的频域GC
% (设工作目录为 prj_neuron_gc/ )
X = gdata(1e5, 3, 5);                  % 计算一个五元三阶的 AR 序列, 长度 len=1e5
wGc = singleGrangerF(X, [1 3], [2]);   % 计算从群组 [1 3] 到 2 的频域 GC
sprintf('GC=%19.16f', mean(wGc))

% 计算两两的频域条件GC
wGC = nGrangerF(X, 10);
g1 = mean(wGC, 3)
g2 = nGrangerT(X, 10)                  % 核对结果
g1-g2
gc_prob_nonzero(abs(g1-g2), 10, 1e5)>0.95   % 相差的量级应与零GC相当

