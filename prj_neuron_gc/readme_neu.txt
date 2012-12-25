这里对工程 nGranger/prj_neuron_gc 里的文件进行一个简单的说明

-------------------------------------------------
变量名约定
simu_time  模拟时间, 单位 ms
stv     采样间隔, 单位 ms
len     数据长度
p       变量数

X       电压(或SpikeTrain等)序列
ISI     平均放电间隔, 单位 ms
ras     放电时刻数据

netstr  网络结构代号 (储存在 prj_neuron_gc/network)
ps	poisson 输入强度
pr	poisson 输入频率(次/ms)
scee	连接强度, 兴奋到兴奋

目录 network 列举了一些可能用到的神经网络的邻接矩阵. gendata_neu.m 中会用到.

-------------------------------------------------

analyse_neuron.m
	示例文件
	
analyse_GC_simple
	示例文件

pic_experiment_all.m
	详细分析一个神经元数据的各个方面.

gendata_neu.m
	[X, ISI, ras] = gendata_neu(netstr, scee, pr, ps, simu_time, stv, extpara)
	神经元数据计算的接口函数. 可作为使用示例.
	extpara 可以填入额外传给计算函数的参数(字符串).
	如果使用 Windows, 注意修改其中的命令行.

pic_experiment_all.m   给出 GC 计算结果的一个较全面的分析报告. (约 195 张图)
pic_experiment_vt.m    显示某一段电压和残差的图像.
pic_experiment_pdf.m   电压和残差的分布
pic_experiment_cov.m   电压和残差的协方差/相关系数
pic_experiment_spike_trigged.m  放电触发的电压和残差


average_volt.m
        Y = average_volt(output_name, id)
        计算群组算术平均. output_name 可以是数据文件路径, 也可以是 X 本身. 向量 id 储存群组编号.
        例如 X 是 9 变量数据(p=9), id = [1,1,1,2,2,2,3,3,3], Y = average_volt(X, id)
        则 Y 是 3 变量数据, 其中 Y(1,:) 是 X(1:3,:) 在每个时刻上的平均.

spikeTriggerAve.m
        [tg_ave, s_rel_time] = spikeTriggerAve(r, s, ras, X, ana_len, stv)
	以放电时刻作为参考点取平均. (所谓的放电触发)

SpikeTrain.m
        [st, fr_cnt] = SpikeTrain(ras, len, neuron_id, bound_ignore, stv, cal_mode)
	生成 Spike Train 序列.

scan_worker/*
        参数扫描

