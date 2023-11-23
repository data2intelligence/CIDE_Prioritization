#!/usr/bin/env python
import os, sys, pathlib, re, pandas, numpy, warnings
import statsmodels.api as sm

from scipy import stats
from glob import glob
from statsmodels.stats.multitest import multipletests
from statsmodels.tools.sm_exceptions import MissingDataError

base_path = pathlib.Path(__file__).parent.absolute()
base_path = os.path.dirname(base_path)
data_path = os.path.join(base_path, 'data')


def matrix_row_stat(f, null_thres, included = None):
    ICB_result = pandas.read_csv(f, sep='\t', index_col=0)
    ICB_result = ICB_result.loc[ICB_result.isnull().mean(axis=1) <= null_thres]
    
    if included is not None:
        ICB_result = ICB_result.loc[ICB_result.index.intersection(included)]
        print('after restriction', ICB_result.shape)
    
    if ICB_result.shape[0] == 0: return None, None
    
    med = ICB_result.median(axis=1)
    med.name = 'med'
    
    p = ICB_result.apply(lambda arr: stats.wilcoxon(arr.loc[~arr.isnull() & (arr != 0)])[1], axis=1)
    p.name = 'p'
    
    result = pandas.concat([med, p], axis=1, join='inner')
    result['FDR'] = multipletests(result['p'], method="fdr_bh")[1]
    
    return ICB_result, result



def prioritize_genes(gene_set, pthres = 0.05, qthres = 0.05, null_thres=0.05, fthres = 2, count_thres = 3, med_thres = 0.5):
    f = os.path.join(data_path, 'merge_immunotherapy.expression.gz')
    
    output = gene_set + '.rank'
    
    included = set()
    fin = open(gene_set)
    for l in fin:
        included.update(l.strip().split('\t'))
    fin.close()
    
    ICB_result, result = matrix_row_stat(f, null_thres, included)
    
    if ICB_result is None:
        print('nothing to rank')
        return
    
    result.to_excel(output + '.stat.xlsx')
    
    flag = (result['p'] < pthres) & (result['FDR'] < qthres)
    result = result.loc[flag]
    
    # with at least one significant hit and no reverse hit
    neg_count = (ICB_result < -fthres).sum(axis=1)
    pos_count = (ICB_result > fthres).sum(axis=1)
    
    flag = (result['med'] < -med_thres) & (neg_count - pos_count >= count_thres)
    result_anti = result.loc[flag].sort_values('med')
    print(result_anti.shape[0], 'negative genes')
    
    flag = (result['med'] > med_thres) & (pos_count - neg_count >= count_thres)
    result_pro = result.loc[flag].sort_values('med', ascending=False)
    print(result_pro.shape[0], 'positive genes')
    
    result = pandas.concat([result_anti, result_pro])
    
    ICB_result = ICB_result.loc[result.index]
    
    writer = pandas.ExcelWriter(output + '.xlsx', engine='xlsxwriter')
    
    result.loc[result['med'] < 0].to_excel(writer, sheet_name='Negative')
    result.loc[result['med'] > 0].to_excel(writer, sheet_name='Positive')
    
    format_number = writer.book.add_format({'num_format': '#,##0.000', 'align': 'center'})
    format_stat = writer.book.add_format({'num_format': '0.00E+00', 'align': 'center'})
    
    for worksheet in writer.sheets:
        worksheet = writer.sheets[worksheet]
        
        worksheet.set_column(1, 1, None, format_number)
        worksheet.set_column(2, 3, None, format_stat)
        
        worksheet.set_zoom(200)
        worksheet.freeze_panes(1, 1)
    
    writer.close()


def main():
    f = sys.argv[1]
    
    if len(sys.argv) < 2 or not os.path.exists(f):
        sys.stderr.write('Please input a gene set\n')
        return 1
    
    prioritize_genes(f)
    
    return 0

if __name__ == '__main__': main()
