GENE_ID_COL_NAME = 'gene_id'
RESULT_TYPE_PVAL = 'pval'
RESULT_TYPE_PROB = 'prob'
# list of (tool, key_col_in_tool_result, tool_result_type)
TOOL_SIG_COL_INFO = [('coloc', 'overall_H4', RESULT_TYPE_PROB),
                     ('fastenloc', 'GLCP', RESULT_TYPE_PROB),
                     ('smr', 'p_SMR', RESULT_TYPE_PVAL),
                     ('predixcan', 'pvalue', RESULT_TYPE_PVAL),
                     ('ecaviar', 'clpp', RESULT_TYPE_PROB),
                     ('twas', 'TWAS.P', RESULT_TYPE_PVAL)
                     ]
AVG_RANKING_COL_NAME = 'avg_ranking'
