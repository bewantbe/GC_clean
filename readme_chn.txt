积分放电模型神经元的格兰杰因果关系(Granger causality, GC)分析程序

安装
    (在命令行下)
    cd /home/user1/matlabtools/      # 假定这是目标安装目录
    git clone https://github.com/bewantbe/GC_clean.git
    cd GC_clean
    matlab -r "install; exit"        # 编译其中的 mex 文件

    至此安装完毕. 如果不编译 mex 文件, 基本的功能仍能使用, 如计算时域GC. 用到的 mex 文件参见 startup.m 的变量 cppfiles.


使用方法
    假定你的 matlab 工作目录是 /home/user1/GC_work/
    新建文件(如果没有的话) /home/user1/GC_work/startup.m
    添加内容:

        run('/home/user1/matlabtools/GC_clean/startup.m');

    其中 /home/user1/matlabtools 是前述的 GC_clean 代码存放目录.
    下次以 /home/user1/GC_work/ 为工作目录启动 matlab 时将会自动添加 GC_clean 的代码目录(addpath)并检查 mex 文件是否最新.
    如果使用的是 Octave, 那么新建 /home/user1/GC_work/.octaverc 并添加同样的内容.

    prj_neuron_gc/analyse_GC_simple.m 是示例文件. 其中用到的可执行文件 raster_tuning (或 raster_tuning.exe), 需要安装 GLUT (例如 freeglut3).
    若 raster_tuning 不适合你的操作系统, 可自行解压 IFsimu_release_2.1.1.zip, 然后编译, 编译好的文件替换 raster_tuning.

    其它详情请看
        GCcal/readme_GCcal.txt
        prj_neuron_gc/readme_neu.txt


常见问题
    若使用 Windows 的记事本查看, 可能有许多文件的换行有问题, 这时请更换你的文本查看器.
    作者主要使用 GNU Octave, 其语法与 Matlab 并不完全相同. 如遇到与 Matlab 不兼容的情况, 请告知 xyy82148在gmail.com. ("在"替换为 @)
    其它关于本代码的建议和Bug也请报告到前述邮箱.


获取最新代码
    hg 版本:
	http://code.google.com/p/granger-causality-neuron/
    git 版本:
	https://github.com/bewantbe/GC_clean.git

    更新记录参见
	http://code.google.com/p/granger-causality-neuron/source/list
	https://github.com/bewantbe/GC_clean/commits/mainline

一般建议总是更新到最新的代码, 已有的代码文件会尽可能保持向后兼容.
