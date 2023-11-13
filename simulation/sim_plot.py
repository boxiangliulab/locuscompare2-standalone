import os

import sim_result_plot


def batch_plot(target_dir, plot_type='ROC'):
    print(f'target dir is {target_dir}')
    if target_dir.startswith('.'):
        raise ValueError(f'The target directory must be absolute path!')
    print(f'started')
    for f in os.listdir(target_dir):
        if not f.startswith('hsq0.'):
            continue
        genetic_model = f.split('_')[1]
        print(f'processing {f}, genetic_model is {genetic_model}')
        if genetic_model == "h20709":
            causal_type = 4
        elif genetic_model == "h20407":
            causal_type = 3
        elif genetic_model == "h2004":
            causal_type = 2
        else:
            causal_type = 0
        generated_file = None
        result_dir = os.path.join(target_dir, f)
        for f1 in os.listdir(result_dir):
            if not f1.startswith('generated_'):
                continue
            generated_file = os.path.join(result_dir, f1)
            break
        if generated_file is None:
            print(f'WARN: There is no generated file in {f}')
            continue
        h1_dir = os.path.join(target_dir, f, 'h1')
        h1_coloc = None
        h1_predixcan = None
        h1_smr = None
        h1_fastenloc = None
        h1_twas = None
        h1_ecaviar = None
        for f2 in os.listdir(h1_dir):
            f2_full_path = os.path.join(h1_dir, f2)
            if f2.startswith('coloc_output_'):
                h1_coloc = f2_full_path
            elif f2.startswith('predixcan_output_'):
                h1_predixcan = f2_full_path
            elif f2.startswith('smr_output_'):
                h1_smr = f2_full_path
            elif f2.startswith('fastenloc_output_'):
                h1_fastenloc = f2_full_path
            elif f2.startswith('twas_output_'):
                h1_twas = f2_full_path
            elif f2.startswith('ecaviar_output_'):
                h1_ecaviar = f2_full_path
        h_gm_dir = os.path.join(target_dir, f, genetic_model)
        hgm_coloc = None
        hgm_predixcan = None
        hgm_smr = None
        hgm_fastenloc = None
        hgm_twas = None
        hgm_ecaviar = None
        for f3 in os.listdir(h_gm_dir):
            f3_full_path = os.path.join(h_gm_dir, f3)
            if f3.startswith('coloc_output_'):
                hgm_coloc = f3_full_path
            elif f3.startswith('predixcan_output_'):
                hgm_predixcan = f3_full_path
            elif f3.startswith('smr_output_'):
                hgm_smr = f3_full_path
            elif f3.startswith('fastenloc_output_'):
                hgm_fastenloc = f3_full_path
            elif f3.startswith('twas_output_'):
                hgm_twas = f3_full_path
            elif f3.startswith('ecaviar_output_'):
                hgm_ecaviar = f3_full_path
        if plot_type == 'ROC':
            roc_path = os.path.join(result_dir, f'{plot_type}_H1_{genetic_model}_OPT.png')
            sim_result_plot.plot_all_against_ensemble_roc(generated_file,
                                                          {'coloc': h1_coloc,
                                                           'smr': h1_smr,
                                                           'fastenloc': h1_fastenloc,
                                                           'predixcan': h1_predixcan,
                                                           'ecaviar': h1_ecaviar,
                                                           'twas': h1_twas},
                                                          {'coloc': hgm_coloc,
                                                           'smr': hgm_smr,
                                                           'fastenloc': hgm_fastenloc,
                                                           'predixcan': hgm_predixcan,
                                                           'ecaviar': hgm_ecaviar,
                                                           'twas': hgm_twas},
                                                          sec_causal_type=causal_type,
                                                          output_figure_path=roc_path,
                                                          output_dir=result_dir)
        elif plot_type == 'PRC':
            prc_path = os.path.join(result_dir, f'{plot_type}_H1_{genetic_model}_OPT.png')
            sim_result_plot.plot_all_against_ensemble_prc(generated_file,
                                                          {'coloc': h1_coloc,
                                                           'smr': h1_smr,
                                                           'fastenloc': h1_fastenloc,
                                                           'predixcan': h1_predixcan,
                                                           'ecaviar': h1_ecaviar,
                                                           'twas': h1_twas},
                                                          {'coloc': hgm_coloc,
                                                           'smr': hgm_smr,
                                                           'fastenloc': hgm_fastenloc,
                                                           'predixcan': hgm_predixcan,
                                                           'ecaviar': hgm_ecaviar,
                                                           'twas': hgm_twas},
                                                          sec_causal_type=causal_type,
                                                          output_figure_path=prc_path,
                                                          output_dir=target_dir)
        elif plot_type == 'BAR':
            bar_path = os.path.join(result_dir, f'{plot_type}_H1_{genetic_model}_OPT.png')
            sim_result_plot.plot_bar(generated_file,
                                     {'coloc': h1_coloc,
                                      'smr': h1_smr,
                                      'fastenloc': h1_fastenloc,
                                      'predixcan': h1_predixcan,
                                      'ecaviar': h1_ecaviar,
                                      'twas': h1_twas},
                                     {'coloc': hgm_coloc,
                                      'smr': hgm_smr,
                                      'fastenloc': hgm_fastenloc,
                                      'predixcan': hgm_predixcan,
                                      'ecaviar': hgm_ecaviar,
                                      'twas': hgm_twas},
                                     sec_causal_type=causal_type,
                                     output_figure_path=bar_path)
        elif plot_type == 'VENN':
            venn_h1_path = os.path.join(result_dir, f'{plot_type}_H1.png')
            venn_gm_path = os.path.join(result_dir, f'{plot_type}_{genetic_model}.png')
            sim_result_plot.plot_venn(generated_file,
                                      {'coloc': h1_coloc,
                                       'smr': h1_smr,
                                       'fastenloc': h1_fastenloc,
                                       'predixcan': h1_predixcan,
                                       'ecaviar': h1_ecaviar,
                                       'twas': h1_twas},
                                      {'coloc': hgm_coloc,
                                       'smr': hgm_smr,
                                       'fastenloc': hgm_fastenloc,
                                       'predixcan': hgm_predixcan,
                                       'ecaviar': hgm_ecaviar,
                                       'twas': hgm_twas},
                                      sec_causal_type=causal_type,
                                      genetic_model=genetic_model,
                                      h1_figure_path=venn_h1_path,
                                      hgm_figure_path=venn_gm_path)
        elif plot_type == 'MEAN':
            mean_std_path = os.path.join(result_dir, f'{plot_type}_H1.png')
            sim_result_plot.plot_mean_sd_combinations_bar(generated_file,
                                                          {'coloc': h1_coloc,
                                                           'smr': h1_smr,
                                                           'fastenloc': h1_fastenloc,
                                                           'predixcan': h1_predixcan,
                                                           'ecaviar': h1_ecaviar,
                                                           'twas': h1_twas},
                                                          {'coloc': hgm_coloc,
                                                           'smr': hgm_smr,
                                                           'fastenloc': hgm_fastenloc,
                                                           'predixcan': hgm_predixcan,
                                                           'ecaviar': hgm_ecaviar,
                                                           'twas': hgm_twas},
                                                          sec_causal_type=causal_type,
                                                          output_figure_path=mean_std_path)
        elif plot_type == 'HEAT':
            try:
                heat_path = os.path.join(result_dir, f'{plot_type}_H1.png')
                sim_result_plot.plot_spearman_heatmap({
                    'coloc': h1_coloc,
                    'smr': h1_smr,
                    'fastenloc': h1_fastenloc,
                    'predixcan': h1_predixcan,
                    'ecaviar': h1_ecaviar,
                    'twas': h1_twas
                }, tested_gene_df=sim_result_plot.retrieve_positive_df(generated_file),
                    output_figure_path=heat_path)
            except Exception as exc:
                print(
                    f'plot corr failed for dir {result_dir}, maybe too many genes are missing in reports of some tools')
        elif plot_type == 'HEAT2':
            try:
                heat_path = os.path.join(result_dir, f'{plot_type}_{genetic_model}.png')
                sim_result_plot.plot_spearman_heatmap({
                    'coloc': hgm_coloc,
                    'smr': hgm_smr,
                    'fastenloc': hgm_fastenloc,
                    'predixcan': hgm_predixcan,
                    'ecaviar': hgm_ecaviar,
                    'twas': hgm_twas
                }, tested_gene_df=sim_result_plot.retrieve_negative_df(generated_file, causal_type),
                    output_figure_path=heat_path)
            except Exception as exc:
                print(
                    f'plot corr failed for dir {result_dir}, maybe too many genes are missing in reports of some tools')
        elif plot_type == 'UPSET':
            mean_std_path = os.path.join(result_dir, f'{plot_type}_H1_{genetic_model}.png')
            sim_result_plot.plot_upset(generated_file,
                                       {'coloc': h1_coloc,
                                        'smr': h1_smr,
                                        'fastenloc': h1_fastenloc,
                                        'predixcan': h1_predixcan,
                                        'ecaviar': h1_ecaviar,
                                        'twas': h1_twas},
                                       {'coloc': hgm_coloc,
                                        'smr': hgm_smr,
                                        'fastenloc': hgm_fastenloc,
                                        'predixcan': hgm_predixcan,
                                        'ecaviar': hgm_ecaviar,
                                        'twas': hgm_twas},
                                       sec_causal_type=causal_type,
                                       output_figure_path=mean_std_path)
        elif plot_type == 'PRF':
            prefix_path = os.path.join(result_dir, f'UNION_{plot_type}_H1_{genetic_model}')
            sim_result_plot.plot_precision_recall_f1_comb_bar(generated_file,
                                                              {'coloc': h1_coloc,
                                                               'smr': h1_smr,
                                                               'fastenloc': h1_fastenloc,
                                                               'predixcan': h1_predixcan,
                                                               'ecaviar': h1_ecaviar,
                                                               'twas': h1_twas},
                                                              {'coloc': hgm_coloc,
                                                               'smr': hgm_smr,
                                                               'fastenloc': hgm_fastenloc,
                                                               'predixcan': hgm_predixcan,
                                                               'ecaviar': hgm_ecaviar,
                                                               'twas': hgm_twas},
                                                              sec_causal_type=causal_type,
                                                              output_figure_prefix=prefix_path,
                                                              typ='UNION')
            prefix_path2 = os.path.join(result_dir, f'INTER_{plot_type}_H1_{genetic_model}')
            sim_result_plot.plot_precision_recall_f1_comb_bar(generated_file,
                                                              {'coloc': h1_coloc,
                                                               'smr': h1_smr,
                                                               'fastenloc': h1_fastenloc,
                                                               'predixcan': h1_predixcan,
                                                               'ecaviar': h1_ecaviar,
                                                               'twas': h1_twas},
                                                              {'coloc': hgm_coloc,
                                                               'smr': hgm_smr,
                                                               'fastenloc': hgm_fastenloc,
                                                               'predixcan': hgm_predixcan,
                                                               'ecaviar': hgm_ecaviar,
                                                               'twas': hgm_twas},
                                                              sec_causal_type=causal_type,
                                                              output_figure_prefix=prefix_path2,
                                                              typ='INTER')
    print(f'completed')


def single_plot(target_dir, genetic_model, plot_type='ROC'):
    print(f'target dir is {target_dir}')
    if target_dir.startswith('.'):
        raise ValueError(f'The target directory must be absolute path!')
    print(f'started')
    if genetic_model == "h20709":
        causal_type = 4
    elif genetic_model == "h20407":
        causal_type = 3
    elif genetic_model == "h2004":
        causal_type = 2
    else:
        causal_type = 0
    generated_file = None
    for f1 in os.listdir(target_dir):
        if not f1.startswith('generated_'):
            continue
        generated_file = os.path.join(target_dir, f1)
        break
    if generated_file is None:
        print(f'WARN: There is no generated file in {target_dir}')
        exit(1)
    h1_dir = os.path.join(target_dir, 'h1')
    h1_coloc = None
    h1_predixcan = None
    h1_smr = None
    h1_fastenloc = None
    h1_twas = None
    h1_ecaviar = None
    for f2 in os.listdir(h1_dir):
        f2_full_path = os.path.join(h1_dir, f2)
        if f2.startswith('coloc_output_'):
            h1_coloc = f2_full_path
        elif f2.startswith('predixcan_output_'):
            h1_predixcan = f2_full_path
        elif f2.startswith('smr_output_'):
            h1_smr = f2_full_path
        elif f2.startswith('fastenloc_output_'):
            h1_fastenloc = f2_full_path
        elif f2.startswith('twas_output_'):
            h1_twas = f2_full_path
        elif f2.startswith('ecaviar_output_'):
            h1_ecaviar = f2_full_path
    h_gm_dir = os.path.join(target_dir, genetic_model)
    hgm_coloc = None
    hgm_predixcan = None
    hgm_smr = None
    hgm_fastenloc = None
    hgm_twas = None
    hgm_ecaviar = None
    for f3 in os.listdir(h_gm_dir):
        f3_full_path = os.path.join(h_gm_dir, f3)
        if f3.startswith('coloc_output_'):
            hgm_coloc = f3_full_path
        elif f3.startswith('predixcan_output_'):
            hgm_predixcan = f3_full_path
        elif f3.startswith('smr_output_'):
            hgm_smr = f3_full_path
        elif f3.startswith('fastenloc_output_'):
            hgm_fastenloc = f3_full_path
        elif f3.startswith('twas_output_'):
            hgm_twas = f3_full_path
        elif f3.startswith('ecaviar_output_'):
            hgm_ecaviar = f3_full_path
    if plot_type == 'ROC':
        roc_path = os.path.join(target_dir, f'{plot_type}_H1_{genetic_model}_OPT.png')
        sim_result_plot.plot_all_against_ensemble_roc(generated_file,
                                                      {'coloc': h1_coloc,
                                                       'smr': h1_smr,
                                                       'fastenloc': h1_fastenloc,
                                                       'predixcan': h1_predixcan,
                                                       'ecaviar': h1_ecaviar,
                                                       'twas': h1_twas},
                                                      {'coloc': hgm_coloc,
                                                       'smr': hgm_smr,
                                                       'fastenloc': hgm_fastenloc,
                                                       'predixcan': hgm_predixcan,
                                                       'ecaviar': hgm_ecaviar,
                                                       'twas': hgm_twas},
                                                      sec_causal_type=causal_type,
                                                      output_figure_path=roc_path,
                                                      output_dir=target_dir)
    elif plot_type == 'PRC':
        prc_path = os.path.join(target_dir, f'{plot_type}_H1_{genetic_model}_OPT.png')
        sim_result_plot.plot_all_against_ensemble_prc(generated_file,
                                                      {'coloc': h1_coloc,
                                                       'smr': h1_smr,
                                                       'fastenloc': h1_fastenloc,
                                                       'predixcan': h1_predixcan,
                                                       'ecaviar': h1_ecaviar,
                                                       'twas': h1_twas},
                                                      {'coloc': hgm_coloc,
                                                       'smr': hgm_smr,
                                                       'fastenloc': hgm_fastenloc,
                                                       'predixcan': hgm_predixcan,
                                                       'ecaviar': hgm_ecaviar,
                                                       'twas': hgm_twas},
                                                      sec_causal_type=causal_type,
                                                      output_figure_path=prc_path,
                                                      output_dir=target_dir)
    elif plot_type == 'BAR':
        bar_path = os.path.join(target_dir, f'{plot_type}_H1_{genetic_model}_OPT.png')
        sim_result_plot.plot_bar(generated_file,
                                 {'coloc': h1_coloc,
                                  'smr': h1_smr,
                                  'fastenloc': h1_fastenloc,
                                  'predixcan': h1_predixcan,
                                  'ecaviar': h1_ecaviar,
                                  'twas': h1_twas},
                                 {'coloc': hgm_coloc,
                                  'smr': hgm_smr,
                                  'fastenloc': hgm_fastenloc,
                                  'predixcan': hgm_predixcan,
                                  'ecaviar': hgm_ecaviar,
                                  'twas': hgm_twas},
                                 sec_causal_type=causal_type,
                                 output_figure_path=bar_path)
    elif plot_type == 'VENN':
        venn_h1_path = os.path.join(target_dir, f'{plot_type}_H1.png')
        venn_gm_path = os.path.join(target_dir, f'{plot_type}_{genetic_model}.png')
        sim_result_plot.plot_venn(generated_file,
                                  {'coloc': h1_coloc,
                                   'smr': h1_smr,
                                   'fastenloc': h1_fastenloc,
                                   'predixcan': h1_predixcan,
                                   'ecaviar': h1_ecaviar,
                                   'twas': h1_twas},
                                  {'coloc': hgm_coloc,
                                   'smr': hgm_smr,
                                   'fastenloc': hgm_fastenloc,
                                   'predixcan': hgm_predixcan,
                                   'ecaviar': hgm_ecaviar,
                                   'twas': hgm_twas},
                                  sec_causal_type=causal_type,
                                  genetic_model=genetic_model,
                                  h1_figure_path=venn_h1_path,
                                  hgm_figure_path=venn_gm_path)
    elif plot_type == 'MEAN':
        mean_std_path = os.path.join(target_dir, f'{plot_type}_H1.png')
        sim_result_plot.plot_mean_sd_combinations_bar(generated_file,
                                                      {'coloc': h1_coloc,
                                                       'smr': h1_smr,
                                                       'fastenloc': h1_fastenloc,
                                                       'predixcan': h1_predixcan,
                                                       'ecaviar': h1_ecaviar,
                                                       'twas': h1_twas},
                                                      {'coloc': hgm_coloc,
                                                       'smr': hgm_smr,
                                                       'fastenloc': hgm_fastenloc,
                                                       'predixcan': hgm_predixcan,
                                                       'ecaviar': hgm_ecaviar,
                                                       'twas': hgm_twas},
                                                      sec_causal_type=causal_type,
                                                      output_figure_path=mean_std_path)
    elif plot_type == 'HEAT':
        try:
            heat_path = os.path.join(target_dir, f'{plot_type}_H1.png')
            sim_result_plot.plot_spearman_heatmap({
                'coloc': h1_coloc,
                'smr': h1_smr,
                'fastenloc': h1_fastenloc,
                'predixcan': h1_predixcan,
                'ecaviar': h1_ecaviar,
                'twas': h1_twas
            }, tested_gene_df=sim_result_plot.retrieve_positive_df(generated_file), output_figure_path=heat_path)
        except Exception as exc:
            print(
                f'plot corr failed for dir {target_dir}, maybe too many genes are missing in reports of some tools')
    elif plot_type == 'HEAT2':
        try:
            heat_path = os.path.join(target_dir, f'{plot_type}_{genetic_model}.png')
            sim_result_plot.plot_spearman_heatmap({
                'coloc': hgm_coloc,
                'smr': hgm_smr,
                'fastenloc': hgm_fastenloc,
                'predixcan': hgm_predixcan,
                'ecaviar': hgm_ecaviar,
                'twas': hgm_twas
            }, tested_gene_df=sim_result_plot.retrieve_negative_df(generated_file, causal_type),
                output_figure_path=heat_path)
        except Exception as exc:
            print(
                f'plot corr failed for dir {target_dir}, maybe too many genes are missing in reports of some tools')
    elif plot_type == 'UPSET':
        mean_std_path = os.path.join(target_dir, f'{plot_type}_H1_{genetic_model}.png')
        sim_result_plot.plot_upset(generated_file,
                                   {'coloc': h1_coloc,
                                    'smr': h1_smr,
                                    'fastenloc': h1_fastenloc,
                                    'predixcan': h1_predixcan,
                                    'ecaviar': h1_ecaviar,
                                    'twas': h1_twas},
                                   {'coloc': hgm_coloc,
                                    'smr': hgm_smr,
                                    'fastenloc': hgm_fastenloc,
                                    'predixcan': hgm_predixcan,
                                    'ecaviar': hgm_ecaviar,
                                    'twas': hgm_twas},
                                   sec_causal_type=causal_type,
                                   output_figure_path=mean_std_path)
    elif plot_type == 'PRF':
        prefix_path = os.path.join(target_dir, f'UNION_{plot_type}_H1_{genetic_model}')
        sim_result_plot.plot_precision_recall_f1_comb_bar(generated_file,
                                                          {'coloc': h1_coloc,
                                                           'smr': h1_smr,
                                                           'fastenloc': h1_fastenloc,
                                                           'predixcan': h1_predixcan,
                                                           'ecaviar': h1_ecaviar,
                                                           'twas': h1_twas},
                                                          {'coloc': hgm_coloc,
                                                           'smr': hgm_smr,
                                                           'fastenloc': hgm_fastenloc,
                                                           'predixcan': hgm_predixcan,
                                                           'ecaviar': hgm_ecaviar,
                                                           'twas': hgm_twas},
                                                          sec_causal_type=causal_type,
                                                          output_figure_prefix=prefix_path,
                                                          typ='UNION')
        prefix_path2 = os.path.join(target_dir, f'INTER_{plot_type}_H1_{genetic_model}')
        sim_result_plot.plot_precision_recall_f1_comb_bar(generated_file,
                                                          {'coloc': h1_coloc,
                                                           'smr': h1_smr,
                                                           'fastenloc': h1_fastenloc,
                                                           'predixcan': h1_predixcan,
                                                           'ecaviar': h1_ecaviar,
                                                           'twas': h1_twas},
                                                          {'coloc': hgm_coloc,
                                                           'smr': hgm_smr,
                                                           'fastenloc': hgm_fastenloc,
                                                           'predixcan': hgm_predixcan,
                                                           'ecaviar': hgm_ecaviar,
                                                           'twas': hgm_twas},
                                                          sec_causal_type=causal_type,
                                                          output_figure_prefix=prefix_path2,
                                                          typ='INTER')
    print(f'completed')


if __name__ == '__main__':
    # if len(sys.argv) < 2:
    #     raise ValueError(f'The target directory is required!')
    # _target_dir = sys.argv[1]

    # batch_plot('/Users/haiyue.meng/Downloads/sim', 'HEAT')
    # batch_plot('/Users/haiyue.meng/Downloads/sim', 'HEAT2')

    single_plot('/Users/haiyue.meng/Downloads/sim/hsq0.2_h2004_rep1', 'h2004', 'HEAT')
