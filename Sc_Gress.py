import scanpy as sc
import numpy as np 
from adjustText import adjust_text
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import shutil
import os

    
######################################################################    
######################################################################
# to be changed by the user, everything else is automatic

dataset_name='GSM4572192'
filepath='GSM4572192_Control1_filtered_feature_bc_matrix.h5'
reference_gene='UMOD'


######################################################################
######################################################################




newdir='figures\\umap' + dataset_name +'\\'
plotdir=dataset_name


if not os.path.exists(newdir):
    os.makedirs(newdir)
    
if not os.path.exists(plotdir):
    os.makedirs(plotdir)




def read_and_prepare_data(filepath, dataset_name,reference_gene):

    dataset=sc.read_10x_h5(filepath)
    dataset.var_names_make_unique
    sc.pl.highest_expr_genes(dataset, n_top=20, )



    dataset.var['mt'] = dataset.var_names.str.startswith('MT-') 
    sc.pp.calculate_qc_metrics(dataset, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)

    sc.pl.violin(dataset, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True)





    adata=dataset
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 10, :]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca(adata, color=reference_gene)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    sc.pl.umap(adata, color=['leiden', reference_gene],use_raw=False,save=dataset_name + '\\Full_data.png')
    adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X


    reference_adata = adata[adata[: , reference_gene].X > 1, :] 
    sc.pp.neighbors(reference_adata, n_neighbors=10, n_pcs=50)
    sc.tl.umap(reference_adata)
    sc.tl.leiden(reference_adata)
    sc.pl.umap(reference_adata, color=['leiden',reference_gene],save=dataset_name + '\\' + reference_gene + '_cells.png')
    sc.tl.rank_genes_groups(reference_adata, 'leiden', method='t-test')
    sc.pl.rank_genes_groups(reference_adata, n_genes=25, sharey=False)
    
    create_peicewise_linear_regressions(reference_adata)





    #specific_cluster=reference_adata[reference_adata.obs['leiden'].isin(['0']),:]
    
def create_peicewise_linear_regressions(reference_adata):


    genes=list(reference_adata.var.index)
    genes

    gene_list=[]
    p_vals=[]
    intercepts=[]
    expression_level=[]
    r_val=[]
    slopes=[]

    for gene in genes:

        reference_subset=reference_adata[: , reference_gene].X
        reference_subset=reference_subset.tolist()


        target_subset=reference_adata[: , gene].X
        target_subset=target_subset.tolist()


        reference_flat= [item for sublist in reference_subset for item in sublist]


        target_flat= [item for sublist in target_subset for item in sublist]

        if len(target_flat)==len(reference_flat)*2:
            target_flat=[1]*len(reference_flat)
        else:
            pass

        exp_level=(reference_adata[reference_adata.obs['leiden'].isin(['0']),:][:,gene].X > 0).mean(0)


        df_pvals = pd.DataFrame({reference_gene:reference_flat,
                          "Gene":target_flat})

        df_pvals=df_pvals[df_pvals['Gene']>1]


        reference_flat=df_pvals[reference_gene].values.tolist()
        target_flat=df_pvals['Gene'].values.tolist()

        try:   
            slope, intercept, r_value, p_value, std_err = stats.linregress(reference_flat,target_flat)
            p_value

        except:
            p_value=1
            intercept=1
            slope=1
            r_value=1

        gene_list.append(gene)
        print(gene)
        p_vals.append(p_value)
        print(p_value)
        intercepts.append(intercept)
        expression_level.append(exp_level)
        r_val.append(r_value)
        slopes.append(slope)

        if p_value < 0.1 and p_value > 0 and exp_level > 0.5:
            sns.scatterplot(reference_flat,target_flat)
            plt.ylabel(gene)
            plt.xlabel(reference_gene)
            plt.title('Correlation of ' + gene + ' with ' + reference_gene + ' p='+ str(round(p_value,6)))

            savetitle=dataset_name + '\\{}_with_{}.png'.format(reference_gene,gene)

            plt.savefig(savetitle)
            plt.clf()

            print(df_pvals)
            #time.sleep(3)


        else:
            pass
        
    create_allgenes_plot(gene_list,p_vals,intercepts,expression_level, r_val, slopes,dataset_name,reference_adata)
        


def create_allgenes_plot(gene_list,p_vals,intercepts,expression_level, r_val, slopes,dataset_name,reference_adata):


    from adjustText import adjust_text


    plt.figure(figsize=(12,8))

    df_regression = pd.DataFrame({"Gene":gene_list,
                          "pval":p_vals,
                                'intercept':intercepts,
                                 'expression':expression_level,
                                 'r_value':r_val,
                                 'slope':slopes})

    df_regression=df_regression.sort_values(by='pval')
    df_regression=df_regression[df_regression['pval']<0.05]
    df_regression=df_regression[df_regression['pval']>0]
    df_regression=df_regression[df_regression['intercept']>0.1]

    exp_vals=df_regression['expression'].values.tolist()

    exp_vals_flat= [item for sublist in exp_vals for item in sublist]

    df_regression['expression']=exp_vals_flat
    df_regression['neglogp']=-np.log(df_regression['pval'])

    df_regression_labs=df_regression[df_regression['expression']>0.5]
    #df_regression_labs=df_regression_labs[df_regression_labs['expression']>0.5]


    df_regression['Positive Correlation'] = df_regression['slope'] > 0
    df_regression['rSq']=df_regression['r_value']**2





    cmap=sns.diverging_palette(250, 30, l=65, center="dark", as_cmap=True)


    sns.scatterplot(x=df_regression['expression'],y=df_regression['neglogp'],hue=df_regression['Positive Correlation'].values.tolist(),s=35)
    plt.legend(title='Correlation Type', loc='upper left', labels=['Positive', 'Negative'])
    plt.ylabel('Linear Regression neglog p-value (AU)')
    plt.xlabel('Proportion of cells also expressing labelled gene')
    plt.ylim(2.5,10)
    #plt.hlines(0.05,0,max(df_regression['expression']),linestyles='--')
    #plt.vlines(0.5,0,max(df_regression['pval']))


    texts=[]
    for x,y,s in zip(df_regression_labs['expression'],df_regression_labs['neglogp'],df_regression_labs['Gene']):
                texts.append(plt.text(x,y,s,size=10))



    adjust_text(texts,autoalign='x'
               )



    plt.savefig(dataset_name+'\Coexpressed_genes_regression.png')
    df_regression_labs.to_csv(plotdir+'\\' + dataset_name+'Coexpressed_genes.csv')
    
    
    x=df_regression_labs.head(3)
    x=x['Gene'].values.tolist()


    sc.pl.umap(reference_adata, color=[reference_gene,x[0],x[1],x[2]],save=dataset_name+'\\Coexpressed_genes.png')
    
    return df_regression_labs

read_and_prepare_data(filepath, dataset_name,reference_gene)


    
source_dir = newdir
target_dir = plotdir
    
file_names = os.listdir(source_dir)
    
for file_name in file_names:
    shutil.move(os.path.join(source_dir, file_name), target_dir)
    
shutil.rmtree(newdir)