import pandas as pd
import scipy.stats
from statistics import mean
import numpy as np
from matplotlib import pyplot as plt
from lifelines import CoxPHFitter
from sklearn.linear_model import LogisticRegression
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import seaborn as sns
import math
import os
from scipy.stats import spearmanr


def scatterplot_all_iter(data, met, mode):
    test = data.loc[data['feature'] == met]
    plt.rcParams['figure.figsize'] = [8, 6]
    if mode == 'iteration':
        test.plot.scatter(x='actual_rank', y='predicted_rank', c='iteration', colormap='twilight')
    if mode == 'censor':
        test.plot.scatter(x='actual_rank', y='predicted_rank', c='censor')
    plt.axline((0, 0), slope=1)
    plt.title(met)
    plt.savefig(f'{plots_dir}/scatterplot_{met}_all_iter_{mode}.pdf')
    plt.close()
    rho, pval = spearmanr(test['actual_rank'], test['predicted_rank'])
    print(f"spearman's rho = {rho}, pval = {pval}")

# quick wrapper for multiple-testing correction
def correct_for_mult(data):
    data['p_adj'] = multipletests(pvals=data['pval'], method="fdr_bh", alpha=0.05)[1]
    return data

def load_data():
    features = set()
    nrows = 0
    samples = []
    batch_names = []
    dirlist = os.listdir(raw_data_dir)
    dirlist.sort()
    for batch_idx, fpath in enumerate(dirlist):
        if 'csv' not in fpath:
            continue
        df = pd.read_csv(f'{raw_data_dir}/{fpath}', header=0, index_col=0).T
        features.update(df.columns)
        samples.extend(df.index)
        batch_names.append(fpath)
        nrows += df.shape[0]
    feature_map = {s: i for i, s in enumerate(features)}
    sample_map = {s: i for i, s in enumerate(samples)}
    batch_map = {s: i for i, s in enumerate(batch_names)}
    data = np.full((nrows, len(features)), np.nan)
    batch_index_vector = np.zeros(nrows, dtype=int)
    sidx = 0
    batch_idx = 0
    for fpath in dirlist:
        if 'csv' not in fpath:
            continue
        df = pd.read_csv(f'{raw_data_dir}/{fpath}', header=0, index_col=0).T
        for feature in df.columns:
            fidx = feature_map[feature]
            data[sidx:sidx + df.shape[0], fidx] = df[feature].values
        batch_index_vector[sidx:sidx + df.shape[0]] = batch_idx
        sidx += df.shape[0]
        batch_idx += 1
    s_test_prop_settings = pd.DataFrame(list(zip(batch_names, proportions)), columns=["batches", "s_test_prop"])
    print("The target dataset is assigned as follows ('1' indicates the target):")
    print(s_test_prop_settings)

    return data, batch_index_vector, feature_map, sample_map, batch_map, s_test_prop_settings


def tic_normalization(data, batch_index_vector):
    normalized_data = np.copy(data)
    n_batches = batch_index_vector.max() + 1
    for bidx in range(n_batches):
        batch_rows = np.arange(data.shape[0])[batch_index_vector == bidx]
        batch = data[batch_rows]
        nan_mask = np.isnan(batch)
        missing = np.all(nan_mask, axis=0)
        min_batch = np.nanmin(batch)
        for row in range(batch.shape[0]):
            n_censored = np.sum(nan_mask[row]) - np.sum(missing)
            row_tic = np.nansum(batch[row]) + 0.5 * n_censored * min_batch
            batch[row, :] = batch[row, :] / row_tic
        batch = batch / np.nansum(batch, axis=1, keepdims=True)
        normalized_data[batch_rows] = batch
    return normalized_data


def convert_to_ranks(data, batch_index_vector):
    n_batches = batch_index_vector.max() + 1
    for bidx in range(n_batches):
        batch_rows = np.arange(data.shape[0])[batch_index_vector == bidx]
        batch = data[batch_rows]
        missing = np.all(np.isnan(batch), axis=0)
        for col in range(data.shape[1]):
            if missing[col]:
                continue
            missing_in_col = np.isnan(batch[:, col])
            rows = batch_rows[~missing_in_col]
            vals = data[rows, col]
            normalizer = len(batch) + 1
            order = vals.argsort()
            data[rows[order], col] = np.arange(missing_in_col.sum() + 1, 1 + len(batch)) / normalizer
            data[batch_rows[missing_in_col], col] = (missing_in_col.sum() + 1) * 0.5 / normalizer
    return data

def get_key(val):
    for key, value in met_map.items():
        if val == value:
            return key
    return "key doesn't exist"

if __name__ == '__main__':
    # Panel A
    ###############
    # Random seed
    seed = 42
    # Input data directory
    raw_data_dir = "MET_raw"
    # Results directory
    results_dir = "results_MET_RNA_imputation/tcga"
    # Whether or not TIC normalization is enabled
    normalize = True
    proportions = [0, 0, 0]
    # Set seed so there is not too much variance in trials
    np.random.seed(seed)
    # Load and preprocess data
    data, batch_index_vector, feature_map, sample_map, batch_map, s_test_prop_settings = load_data()
    features_from_index = {v: k for k, v in feature_map.items()}
    feature_names = [features_from_index[i] for i in range(len(feature_map))]
    samples_from_index = {v: k for k, v in sample_map.items()}
    sample_names = [samples_from_index[i] for i in range(len(sample_map))]
    if normalize:
        data = tic_normalization(data, batch_index_vector)
    ranked_data = convert_to_ranks(data, batch_index_vector)
    met_ranked_data = pd.DataFrame(ranked_data, index=sample_names, columns=feature_names)
    met_ranked_data.to_csv(f'{results_dir}/met_ranked_data.csv')
    # Load the confident predicted metabolites
    metabolites = pd.read_csv(f'results_MET_RNA/reproducibly_well_predicted_met_3.csv', header=0, index_col=0)
    len_met = metabolites.shape[0]
    met = []
    met.extend(metabolites.index)
    met_map = {s: f'{i}a' for i, s in enumerate(met)}
    # TCGA and RC12
    sub_dir = 'tcga'
    results_dir = f'results_MET_RNA_imputation/{sub_dir}'
    clinical_rc12_path = f'{results_dir}/rc12_staginginfo.xlsx'
    clinical_tcga_path = f'{results_dir}/KIRC.clin.merged.picked.txt'
    sample_type_path = f'{results_dir}/sample_type.csv'
    len_samples_tcga = 606
    len_samples_rc12 = 280  # all rc12 samples

    # Load the ranked MET_RNA matrix
    met_ranked_data = pd.read_csv(f'{results_dir}/met_ranked_data.csv', header=0, index_col=0)
    rc12 = met_ranked_data.iloc[:len_samples_rc12, :]
    rc12 = rc12.dropna(axis=1, how='all')
    clinical_rc12 = pd.read_excel(clinical_rc12_path, index_col="SAMPLE NAME")
    rc12 = pd.concat([rc12,clinical_rc12],axis = 1)

    # Load the predicted MET_RNA matrix
    ranked_predictions = pd.read_csv(f'{results_dir}/predicted_data_ave.csv', header=0, index_col=0)
    tcga = ranked_predictions.iloc[:len_samples_tcga, :1573]
    sample_type_tcga = pd.read_csv(sample_type_path, header=0, index_col=0)
    tcga = pd.concat([tcga,sample_type_tcga],axis = 1)

    # Load the clinical data
    clinical_tcga = df = pd.read_csv(clinical_tcga_path, delimiter="\t",header=0, index_col=0).T
    ################# Wilcoxon rank sum test (= Mann Whitney U Test ) ----- Tumor vs Normal Samples
    # rc12
    test = rc12
    met_wilcox_df = pd.DataFrame({'log2fc': [], 'p': []})
    for i in list(test.iloc[:, :877].columns):
        tumor = test.loc[test['TISSUE TYPE'] == 'TUMOR', i]
        normal = test.loc[(test['TISSUE TYPE'] == 'NORMAL', i)]
        met_wilcox_df.loc[i, 'log2fc'] = tumor.mean() - (normal.mean())
        met_wilcox_df.loc[i, 'p'] = scipy.stats.ranksums(tumor, normal).pvalue

    # tcga
    test = tcga
    met_wilcox_df_tcga = pd.DataFrame({'log2fc_tcga': [], 'p_tcga': []})
    for i in list(test.iloc[:, :1573].columns):
        tumor = test.loc[test['TN'] == 'Tumor', i]
        normal = test.loc[(test['TN'] == 'Normal', i)]
        met_wilcox_df_tcga.loc[i, 'log2fc_tcga'] = tumor.mean() - (normal.mean())
        met_wilcox_df_tcga.loc[i, 'p_tcga'] = scipy.stats.ranksums(tumor, normal).pvalue
    a = set()
    a.update(met_wilcox_df_tcga.index)
    b = set()
    b.update(met_wilcox_df.index)
    common_met = a.intersection(b)
    met_wilcox_df = met_wilcox_df.loc[common_met]
    met_wilcox_df_tcga = met_wilcox_df_tcga.loc[common_met]

    # multipletest correction
    met_wilcox_df['p_adj'] = multipletests(pvals=met_wilcox_df['p'], method="fdr_bh", alpha=0.05)[1]
    met_wilcox_df_tcga['p_adj_tcga'] = multipletests(pvals=met_wilcox_df_tcga['p_tcga'], method="fdr_bh", alpha=0.05)[1]
    wilcox_tn = pd.concat([met_wilcox_df, met_wilcox_df_tcga], axis=1)
    wilcox_tn.to_csv(f'{results_dir}/wilcox_tn_rc12_tcga.csv')

    # Only look at the 262 confidently predicted metabolites
    wilcox_tn = pd.read_csv(f'{results_dir}/wilcox_tn_rc12_tcga.csv', header=0, index_col=0)
    a = set()
    a.update(wilcox_tn.index)
    b = set()
    b.update(metabolites.index)
    common_met = a.intersection(b)
    wilcox_tn = wilcox_tn.loc[common_met]
    wilcox_tn['p_adj'] = multipletests(pvals=wilcox_tn['p'], method="fdr_bh", alpha=0.1)[1]
    wilcox_tn['p_adj_tcga'] = multipletests(pvals=wilcox_tn['p_tcga'], method="fdr_bh", alpha=0.1)[1]
    wilcox_tn.to_csv(f'{results_dir}/wilcox_tn_rc12_tcga_110met_mean_difference.csv')

    ############################## Wilcoxon rank sum test (= Mann Whitney U Test ) ----- High vs Low Stages
    # rc12
    test = rc12.loc[rc12['TISSUE TYPE'] == 'TUMOR']
    met_wilcox_df = pd.DataFrame({'log2fc': [], 'p': []})
    for i in list(test.iloc[:, :877].columns):
        low = test.loc[test['PATHOLOGY STAGE - 3'] == 'STAGES I + II', i]
        high = test.loc[(test['PATHOLOGY STAGE - 3'] == 'STAGES III + IV', i)]
        met_wilcox_df.loc[i, 'log2fc'] = high.mean() - (low.mean())
        met_wilcox_df.loc[i, 'p'] = scipy.stats.ranksums(high, low).pvalue

    # tcga
    test = tcga.loc[tcga['TN'] == 'Tumor']
    test['short_brc'] = [brc[0:12].lower() for brc in test.index]
    test = test.drop(index='TCGA-DV-A4W0-05A-11R-A266-07')  # drop the 05A duplicate samples for one patient
    test = pd.merge(test, clinical_tcga, left_on='short_brc', right_on=clinical_tcga.index, validate='one_to_one')
    test.loc[(test['pathologic_stage'] == 'stage i') | (test['pathologic_stage'] == 'stage ii'), 'stage'] = 'STAGES I + II'
    test.loc[
    (test['pathologic_stage'] == 'stage iii') | (test['pathologic_stage'] == 'stage iv'), 'stage'] = 'STAGES III + IV'

    met_wilcox_df_tcga = pd.DataFrame({'log2fc_tcga': [], 'p_tcga': []})
    for i in list(test.iloc[:, :1573].columns):
        low = test.loc[test['stage'] == 'STAGES I + II', i]
        high = test.loc[(test['stage'] == 'STAGES III + IV', i)]
        met_wilcox_df_tcga.loc[i, 'log2fc_tcga'] = high.mean() - (low.mean())
        met_wilcox_df_tcga.loc[i, 'p_tcga'] = scipy.stats.ranksums(high, low).pvalue

    a = set()
    a.update(met_wilcox_df_tcga.index)
    b = set()
    b.update(met_wilcox_df.index)
    common_met = a.intersection(b)
    met_wilcox_df = met_wilcox_df.loc[common_met]
    met_wilcox_df_tcga = met_wilcox_df_tcga.loc[common_met]

    # multipletest correction
    met_wilcox_df['p_adj'] = multipletests(pvals=met_wilcox_df['p'], method="fdr_bh", alpha=0.05)[1]
    met_wilcox_df_tcga['p_adj_tcga'] = multipletests(pvals=met_wilcox_df_tcga['p_tcga'], method="fdr_bh", alpha=0.05)[1]
    wilcox_stage = pd.concat([met_wilcox_df, met_wilcox_df_tcga], axis=1)
    wilcox_stage.to_csv(f'{results_dir}/wilcox_stage_rc12_tcga.csv')

    # Only look at the 262 confidently predicted metabolites
    wilcox_stage = pd.read_csv(f'{results_dir}/wilcox_stage_rc12_tcga.csv', header=0, index_col=0)
    a = set()
    a.update(wilcox_stage.index)
    b = set()
    b.update(metabolites.index)
    common_met = a.intersection(b)
    wilcox_stage = wilcox_stage.loc[common_met]
    wilcox_stage['p_adj'] = multipletests(pvals=wilcox_stage['p'], method="fdr_bh", alpha=0.1)[1]
    wilcox_stage['p_adj_tcga'] = multipletests(pvals=wilcox_stage['p_tcga'], method="fdr_bh", alpha=0.1)[1]
    wilcox_stage.to_csv(f'{results_dir}/wilcox_stage_rc12_tcga_110met_mean_difference.csv')

    # Panel B
    plots_dir = f'{results_dir}/plots'
    actual_pred_res_df = pd.read_csv(f'{results_dir}/actual_vs_predicted_ranks.csv', header=0, index_col=0)
    # Calculate rho and p-values
    cor = actual_pred_res_df.groupby(["f_test_prop", "iteration", "feature"])[['actual_rank', 'predicted_rank']].corr(
        method='spearman').iloc[0::2, -1].to_frame()
    pval = actual_pred_res_df.groupby(["f_test_prop", "iteration", "feature"])[['actual_rank', 'predicted_rank']].corr(
        method=lambda x, y: spearmanr(x, y)[1]).iloc[0::2, -1].to_frame()

    # obtain rho for each feature in each iteration
    by_iter_rho = (pd.merge(cor.rename(columns={"predicted_rank": "rho"}),
                            pval.rename(columns={"predicted_rank": "pval"}),
                            how='left',
                            on=['f_test_prop', 'iteration', 'feature'])
                   .groupby(["f_test_prop", "iteration"], group_keys=False)
                   .apply(correct_for_mult)
                   # bin into significant/ns
                   .assign(sig=lambda dataframe: dataframe['p_adj'].map(lambda p_adj: True if p_adj < 0.05 else False))
                   )

    # calculate median rho for each feature
    # once again, z-transform and group appropriately
    g_median_by_feature = (by_iter_rho
                           .assign(z_score=lambda dataframe: dataframe['rho'].map(lambda rho: np.arctanh(rho)))
                           .groupby(["f_test_prop", "feature"], group_keys=False))
    median_rho_feature = (pd.DataFrame({"median_z": g_median_by_feature.median()['z_score'],
                                        "sig_in": g_median_by_feature.sum()['sig'] / g_median_by_feature.size()})
                          .assign(median_rho=lambda dataframe: dataframe['median_z'].map(lambda z: np.tanh(z)))
                          # flag well-predicted metabolites
                          .assign(sig_in_most=(lambda dataframe: dataframe['sig_in']
                                               # note that if a flag for "well-predicted metabolite" is desired, a condition for positive rho must be added here
                                               .map(lambda sig_in: True if sig_in > 0.90 else False)))
                          .reset_index()
                          )
    # Bar plot
    bar_df = (median_rho_feature[median_rho_feature["f_test_prop"] == min(median_rho_feature["f_test_prop"])]
              .sort_values(by=["median_rho"], ascending=False)
              .assign(color=lambda df: df['sig_in_most'].map(lambda wp: 'red' if wp else 'grey'))
              )

    plt.rcParams['figure.figsize'] = [30, 40]
    plt.bar(bar_df["feature"], bar_df["median_rho"], color=bar_df["color"])
    plt.xlabel("metabolite")
    plt.ylabel("median rho")
    plt.title("median rho for metabolites at X test proportion")
    plt.xticks(rotation=-90)
    plt.margins(x=0.01)
    plt.savefig(f'{plots_dir}/Bar plot of median rho by metabolites.pdf')
    plt.close()

    # Panel C scatter plot
    # scatter plot for the most well-predicted metabolite
    scatterplot_all_iter(actual_pred_res_df,
                         median_rho_feature.loc[median_rho_feature['median_rho'] == median_rho_feature['median_rho'].
                         max(), 'feature'].values[0], mode='iteration')

    # Panel E
    # Kaplan-Meier plot for 1-methylimidazole acetate
    # top
    # Javelin 101
    sub_dir = 'javelin_101'
    results_dir = f'results_MET_RNA_imputation/{sub_dir}'
    clinical_path = f'Other_RNA/clinical data/{sub_dir}'
    len_samples = 726
    # Load the predicted MET_RNA matrix
    ranked_predictions = pd.read_csv(f'{results_dir}/predicted_data_ave.csv', header=0, index_col=0)
    met_predicted = ranked_predictions.iloc[:len_samples, :]
    met_predicted = met_predicted.iloc[:, :1573]
    met_predicted = met_predicted.loc[:, metabolites.index]
    met_predicted = met_predicted.rename(columns=met_map)
    # Load the clinical data
    clinical = pd.read_excel(f'{clinical_path}.xlsx', sheet_name="all_clinical_data", index_col='ID')
    # Merge the data
    data = pd.concat([met_predicted, clinical], axis=1)
    data = data.iloc[:len_samples, :]
    data = data.rename({'PFS_P': 'PFS_MO'}, axis=1)
    data.loc[data['PFS_P_CNSR'] == 0, 'PFS_EVENT'] = 1
    data.loc[data['PFS_P_CNSR'] == 1, 'PFS_EVENT'] = 0
    cols = []
    cols.extend(met_map.values())
    cols.extend(['AGE', 'SEX', 'PFS_MO', 'PFS_EVENT', 'TRT01P'])
    survival_df_add = data[cols]
    survival_df_add['Dataset'] = 'javelin_101'
    survival_df_add_C = survival_df_add.loc[survival_df_add['TRT01P'] == 'Sunitinib']
    print(survival_df_add_C['26a'].median())  # 0.4711141678129299

    test = survival_df_add_C
    T1 = test[test['26a'] < survival_df_add_C['26a'].median()]['PFS_MO']
    T2 = test[test['26a'] >= survival_df_add_C['26a'].median()]['PFS_MO']
    E1 = test[test['26a'] < survival_df_add_C['26a'].median()]['PFS_EVENT']
    E2 = test[test['26a'] >= survival_df_add_C['26a'].median()]['PFS_EVENT']
    kmf = KaplanMeierFitter(label="High 1-methylimidazoleacetate")
    kmf.fit(T2, E2)
    kmf.plot(show_censors=True, ci_show=False)
    print(kmf.median_survival_time_)
    kmf = KaplanMeierFitter(label="Low 1-methylimidazoleacetate")
    kmf.fit(T1, E1)
    kmf.plot(show_censors=True, ci_show=False)
    print(kmf.median_survival_time_)

    results = logrank_test(T1, T2, event_observed_A=E1, event_observed_B=E2)
    print(results.p_value)
    plt.xlabel("Progression Free Survival(months)")
    plt.ylabel("Survival Probability")
    plt.show()

    plt.savefig('results_MET_RNA_imputation/tcga/plots/KM_1-methylimidazoleacetate_Javelin101.pdf')

    # bottom
    # Comparz
    sub_dir = 'Comparz'
    results_dir = f'results_MET_RNA_imputation/{sub_dir}'
    clinical_path = f'Other_RNA/clinical data/{sub_dir}'
    len_samples = 412
    # Load the predicted MET_RNA matrix
    ranked_predictions = pd.read_csv(f'{results_dir}/predicted_data_ave.csv', header=0, index_col=0)
    met_predicted = ranked_predictions.iloc[:len_samples, :]
    met_predicted = met_predicted.iloc[:, :1573]
    met_predicted = met_predicted.loc[:, metabolites.index]
    met_predicted = met_predicted.rename(columns=met_map)
    # Load the clinical data
    clinical = pd.read_excel(f'{clinical_path}.xlsx', sheet_name="Comparz.FOLH1", index_col="RNASampleID")
    clinical = clinical.iloc[25:, :]
    # Merge the data
    data = pd.concat([met_predicted, clinical], axis=1)
    data = data.iloc[:len_samples, :]
    data = data[data['SRVCFLCD'].notna()]
    data = data[data['SRVMO'].notna()]
    data = data.rename({'SRVCFLCD': 'OS_EVENT'}, axis=1)
    data = data.rename({'SRVMO': 'OS_MO'}, axis=1)
    data = data.rename({'PFSCFLCD': 'PFS_EVENT'}, axis=1)
    data = data.rename({'PFSMO': 'PFS_MO'}, axis=1)
    cols = []
    cols.extend(met_map.values())
    cols.extend(['AGE', 'SEX', 'PFS_MO', 'PFS_EVENT', 'OS_EVENT', 'OS_MO', 'TRTGRP'])
    survival_df_add = data[cols]
    survival_df_add['Dataset'] = 'Comparz'
    survival_df_add_C = survival_df_add.loc[survival_df_add['TRTGRP'] == 'sunitinib']
    print(survival_df_add_C['26a'].median())  # 0.5211864406779662

    test = survival_df_add_C
    T1 = test[test['26a']<survival_df_add_C['26a'].median()]['PFS_MO']
    T2 = test[test['26a']>=survival_df_add_C['26a'].median()]['PFS_MO']
    E1 = test[test['26a']<survival_df_add_C['26a'].median()]['PFS_EVENT']
    E2 = test[test['26a']>= survival_df_add_C['26a'].median()]['PFS_EVENT']
    kmf = KaplanMeierFitter(label="High 1-methylimidazoleacetate")
    kmf.fit(T2, E2)
    kmf.plot(show_censors=True,ci_show=False)
    print(kmf.median_survival_time_)
    kmf = KaplanMeierFitter(label="Low 1-methylimidazoleacetate")
    kmf.fit(T1, E1)
    kmf.plot(show_censors=True,ci_show=False)
    print(kmf.median_survival_time_)

    results = logrank_test(T1, T2, event_observed_A=E1, event_observed_B=E2)
    print(results.p_value)
    plt.xlabel("Progression Free Survival(months)")
    plt.ylabel("Survival Probability")
    plt.show()

    plt.savefig('results_MET_RNA_imputation/tcga/plots/KM_1-methylimidazoleacetate_COMPARZ.pdf')

