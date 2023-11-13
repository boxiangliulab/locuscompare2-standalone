# Astle heatmap
相关处理在simulation/real_plot.py中

1.  前置处理(可跳过, 结果已经在下载文件中了)
```
gwas_tigssue
calc_overview_data
```

2. plot
```commandline
plot_clustermap_with_errorbar
plot_clustermap_with_ci
```

# sem sim data ROC curve
相关处理在simulation/sem_plot.py中

1. 前置处理(可跳过, 结果已经在下载文件中了)
```commandline
merge_result
```

2. plot
```commandline
plot_all_against_ensemble_roc
```

# precision recall bar plot

相关处理在simulation/sem_plot.py中.
* 下载文件的sem_sim_data.zip中有2个文件夹:split_no_limit和original_smr_no_p_limit.
* split_no_limit是把所有的threshold都去掉后执行colocalization的结果.
* original_smr_no_p_limit中只包含原始smr(把threshold去掉)的结果.
* split_no_limit/smr_merged.tsv就是original_smr_no_p_limit中smr结果的合并.
* split_no_limit/ENSG*/smr_output_*.tsv.gz这里面是自己实现的single snp smr在每个基因里的结果(现在plot没有用这个).


1. 前置处理(可跳过, 结果已经在下载文件中了)
```commandline
merge_result
```

2. 不同个数工具precision recall union/intersection bar plot
```commandline
plot_precision_recall_f1_comb_bar
```

3. 每个工具precision/recall和majority vote precision recall比对
```commandline
plot_tools_compare_precision_recall_f1_bar
```
如果majority vote需要改m值, 就把下面代码中的3换成其他值
```commandline
merged_df['mvote_positive'] = merged_df[mark_cols].sum(axis=1) >= 3
```
