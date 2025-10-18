
现在帮我完成这个可视化软件，我已经完成了大部分需要的函数，现在需要你做以下事情：

1. 将app.R分成ui.R和server.R
2. 使用我写的函数库functions.R完成server.R。以下是我对这个软件功能的描述，软件要根据用户不同的输入文件以及选项调整不同的输出项目：
   1. 用户输入LMM的rds结果文件（药物反应数据，lineage barcode），scrna,scatac,。对于非multiomic数据有rna seq的或者atac数据的lineage barcode与cell barcode mapping的表，multiomic数据统一有一张mapping表。程序需要根据不同的输入调整输出，例如用户只上传了scrna或者scatac数据，或者没有上传rna 或atac的mapping表，程序都需要舍弃一些输出。对于LMM的rds结果文件，用户可以选择对哪些数据denoise，并储存在dslt里。同时还要实现用户上传一个另外的药物表达矩阵，选择是否denoise，并在之后的绘图界面可选。
   2. 用户选择好文件之后，并选择需要筛选的barcode后缀，点击一个按钮，程序读取所选的文件。
   3. 用户切换到qc setting，切换到qc setting页面之后，设置好参数之后，输出qc图，程序按照参数subset。然后使用函数工具将药物，rna，atac数据按照mapping表映射到lineage或者sc level，之后对这些数据找cluster以及降维，储存在dslt里。另外对药物反应数据查找archetype。
   4. 之后根据用户选择的cluster方法（louvain或者kmeans）以及降维方法（pca，tsne，umap）还有维度绘制点图。lineage，sc两个level分别绘制rna，atac，以及选择的药物反应数据，三张降维图，再特定药物绘制点热图。然后选择根据cluster，哪种cluster方法，或者archetype来绘制气泡图。还有根据输入的基因以及区域绘制小提琴组图以及feature plot
   5. 再输出一些有用的东西
3. 你需要根据这些描述，根据里面描述的逻辑，组合functions.R里的函数，实现上面的功能。
4. 若现有的代码无法完成描述的功能，你可以添加与修改，但是代码必须简洁，只保留必要代码。
5. 同时，尽量不要改动现有代码的逻辑。