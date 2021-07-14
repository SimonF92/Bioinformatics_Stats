import scanpy as sc
import numpy as np 
from adjustText import adjust_text
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats



GSM457=sc.read_10x_h5('GSM4572193/GSM4572193_Control2_filtered_feature_bc_matrix.h5')

GSM457.var_names_make_unique

sc.pl.highest_expr_genes(GSM457, n_top=20, )



GSM457.var['mt'] = GSM457.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(GSM457, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)

sc.pl.violin(GSM457, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)


reference_gene='UMOD'


adata=GSM457

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
sc.pl.umap(adata, color=['leiden', reference_gene,'SLC12A1'],use_raw=False,save='GSM4572193\\Full_data.png')

adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
umod_adata = adata[adata[: , reference_gene].X > 1, :] 
sc.pp.neighbors(umod_adata, n_neighbors=10, n_pcs=50)


sc.tl.umap(umod_adata)
sc.tl.leiden(umod_adata)

sc.pl.umap(umod_adata, color=['leiden',reference_gene,'SLC12A1','CASR'],save='GSM4572193\\UMOD_cells.png')


sc.tl.rank_genes_groups(umod_adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(umod_adata, n_genes=25, sharey=False)





umod_cluster=umod_adata[umod_adata.obs['leiden'].isin(['0']),:]

genes=list(umod_cluster.var.index)
genes

gene_list=[]
p_vals=[]
intercepts=[]
expression_level=[]
r_val=[]
slopes=[]

for gene in genes:
    
    umod=umod_cluster[: , reference_gene].X
    umod=umod.tolist()


    gene_values=umod_cluster[: , gene].X
    gene_values=gene_values.tolist()


    umod_flat= [item for sublist in umod for item in sublist]


    gene_flat= [item for sublist in gene_values for item in sublist]
    
    if len(gene_flat)==len(umod_flat)*2:
        gene_flat=[1]*len(umod_flat)
    else:
        pass
    
    exp_level=(umod_adata[umod_adata.obs['leiden'].isin(['0']),:][:,gene].X > 0).mean(0)


    df_pvals = pd.DataFrame({reference_gene:umod_flat,
                      "Gene":gene_flat})

    df_pvals=df_pvals[df_pvals['Gene']>1]


    umod_flat=df_pvals[reference_gene].values.tolist()
    gene_flat=df_pvals['Gene'].values.tolist()
    
    try:   
        slope, intercept, r_value, p_value, std_err = stats.linregress(umod_flat,gene_flat)
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
    
    
    
    



plt.figure(figsize=(12,8))

df_regression = pd.DataFrame({"Gene":gene_list,
                      "pval":p_vals,
                            'intercept':intercepts,
                             'expression':expression_level,
                             'r_value':r_val,
                             'slope':slopes})

df_regression=df_regression.sort_values(by='pval')
df_regression=df_regression[df_regression['pval']<0.2]
df_regression=df_regression[df_regression['pval']>0]
df_regression=df_regression[df_regression['intercept']>0.1]

exp_vals=df_regression['expression'].values.tolist()

exp_vals_flat= [item for sublist in exp_vals for item in sublist]

df_regression['expression']=exp_vals_flat


df_regression_labs=df_regression[df_regression['expression']>0.5]
#df_regression_labs=df_regression_labs[df_regression_labs['expression']>0.5]


df_regression['Positive Correlation'] = df_regression['slope'] > 0
df_regression['rSq']=df_regression['r_value']**2


texts=[]
for x,y,s in zip(df_regression_labs['expression'],df_regression_labs['pval'],df_regression_labs['Gene']):
            texts.append(plt.text(x,y,s,size=10))
        

        
adjust_text(texts,precision=0.5,
        expand_text=(1.01, 1.05), expand_points=(1.01, 1.05),
        force_text=(0.01, 0.25), force_points=(0.01, 0.25)
        )


cmap=sns.diverging_palette(250, 30, l=65, center="dark", as_cmap=True)


sns.scatterplot(x=df_regression['expression'],y=df_regression['pval'],hue=df_regression['Positive Correlation'].values.tolist(),s=35)
plt.legend(title='Correlation Type', loc='upper right', labels=['Positive', 'Negative'])
plt.ylabel('Linear Regression p-value (AU)')
plt.xlabel('Relative Expression in Cluster to ' + reference_gene)
#plt.hlines(0.05,0,max(df_regression['expression']),linestyles='--')
#plt.vlines(0.5,0,max(df_regression['pval']))

plt.savefig('GSM4572193\Coexpressed_genes.png')
df_regression_labs.to_csv('GSM4572193\Coexpressed_genes.csv')


x=df_regression_labs.head(3)
x=x['Gene'].values.tolist()


sc.pl.umap(umod_adata, color=[reference_gene,x[0],x[1],x[2]],save='GSM4572193\\Coexpressed_genes.png')