#author::: Simon Fisher s.fisher.1@research.gla.ac.uk


#may have to pip install lqrt

import pandas as pd
import numpy as np
import math
import collections
import scipy.stats
import scipy.stats.distributions as distributions
import numpy as np
import lqrt



df=pd.read_table("orderedRawCounts.tsv")

group1=['GX1403_T0','GX1404_T0','GX1405_T0','GX1406_T0','GX1410_T0','GX1413_T0','GX1414_T0']
group2=['GX1403_T3','GX1404_T3','GX1405_T3','GX1406_T3','GX1410_T3','GX1413_T3','GX1414_T3']
genelist=df.index.tolist()

dfgroup1=None
dfgroup1=pd.DataFrame()
dfgroup1['ID']=genelist




###########uses list comprehension not pandas to track errors################

for sample in group1:
    x=df[sample].values.tolist()  
    
    
    prepped=[] 
    #log+1 counts
    for value in x:
        value=np.log(value+1)
        prepped.append(value)
    
    #calculate sample total reads (log+1 version)
    total=sum(prepped) 
    
    normalised=[]
    for value in prepped:
        #normalise to total
        normvalue= value/total
        normalised.append(normvalue)
        

    
    dfgroup1[sample]=normalised


dfgroup1=dfgroup1.set_index('ID')       
dfgroup1['group1_mean'] = dfgroup1.mean(axis=1)   




#repeat for group 2
dfgroup2=None
dfgroup2=pd.DataFrame()
dfgroup2['ID']=genelist

for sample in group2:
    x=df[sample].values.tolist()  
    
    
    prepped=[] 
    for value in x:
        value=np.log(value+1)
        prepped.append(value)
        
    total=sum(prepped) 
    
    normalised=[]
    for value in prepped:
        normvalue= value/total
        normalised.append(normvalue)
        

    
    dfgroup2[sample]=normalised
    


dfgroup2=dfgroup2.set_index('ID')       
dfgroup2['group2_mean'] = dfgroup2.mean(axis=1)  

    
    
##########################################################################

    
    
#combine into dataframe

normaliseddf= pd.concat([dfgroup1,dfgroup2], axis=1)
#kill rows with zero values
normaliseddf = normaliseddf.replace(0, np.nan)
normaliseddf=normaliseddf.dropna()


group1_mean=  normaliseddf['group1_mean'].values.tolist()    
group1_mean

group2_mean=  normaliseddf['group2_mean'].values.tolist()    
group2_mean

l2fcs=[]

for x, y in zip(group1_mean, group2_mean):
    #does not log transform here because log trans was already performed (this is a bit im not sure on)
    l2fc= (y/x)
    l2fcs.append(l2fc)
    
normaliseddf['l2fc']= l2fcs

normaliseddf=normaliseddf.reset_index()



#performs both students and welchs ttest, treats data as simple array, no correction for multiple testing due to low gene counts to begin with   
normaliseddf['standard t-test']=scipy.stats.ttest_ind(normaliseddf.iloc[:, 1:8], normaliseddf.iloc[:, 9:16], axis=1)[1]
#welchs correction due to possible unequal variance
normaliseddf['welch correction']=scipy.stats.ttest_ind(normaliseddf.iloc[:, 1:8], normaliseddf.iloc[:, 9:16], axis=1,equal_var=False)[1]





##########################################################################
##########################################################################
##########################################################################
# Whole next block calculates the Wald stat, the same test that Deseq2 uses

def WaldTest(x,y):
    
    sampleSize=len(group1)

    X=x
    Y=y

    # Compute their average
    avgX = sum(X)/sampleSize
    avgY = sum(Y)/sampleSize
    


    # Partial steps to compute estimators of linear regression parameters.
    XDiff = [X_i - avgX for X_i in X]
    XDiffSquared = [i*i for i in XDiff]
    YDiff = [Y_i - avgY for Y_i in Y]

    # B1 is the estimator of slope.
    # B0 is the estimator of intercept.
    # r is the estimator of Y given X.
    B1 = sum(x * y for x, y in zip(XDiff, YDiff)) / sum(XDiffSquared)
    B0 = avgY - B1*avgX
    r = lambda x: B0 + B1*x

    # Partial steps to compute Wald Statistic.
    errs = [y - r(x) for x, y in zip(X, Y)]
    errStd = math.sqrt((1/(sampleSize-2))*(sum([err**2 for err in errs])))
    XStd = math.sqrt((1/(sampleSize))*sum([diff**2 for diff in XDiff]))
    stdB1 = errStd / (XStd * math.sqrt(sampleSize))

    # Wald Statistic.
    W = (B1 - 0)/stdB1

    # pvalue of Wald Test of B1 = 0.
    pvalueWald = 2*scipy.stats.norm.cdf(-abs(W))
    return pvalueWald
    x=None
    y=None


    
array= normaliseddf.iloc[:, 1:8].values.tolist()
flat_array = [item for sublist in array for item in sublist]
flat_array



waldvalues=[]
length=len(normaliseddf)
length
i=0

while i <= (length-1):
    
    subset=normaliseddf.loc[i]
    waldvalue = WaldTest(list(subset[1:8]), list(subset[9:16]))
    waldvalues.append(waldvalue)
    i+=1
    
    
##########################################################################
##########################################################################
##########################################################################


mannwhitneys=[]
length=len(normaliseddf)
length
i=0

while i <= (length-1):
    
    subset=normaliseddf.loc[i]
    x= ((scipy.stats.mannwhitneyu(list(subset[1:8]), list(subset[9:16]), use_continuity=False,alternative='two-sided'))[1])
    
    mannwhitneys.append(x)
    i+=1
    

    
################the following test will be slow with low computer specs#####################
'''
lqrts=[]
length=len(normaliseddf)
length
i=0

while i <= (length-1):
    
    subset=normaliseddf.loc[i]
    x=lqrt.lqrtest_ind(list(subset[1:8]), list(subset[9:16]), equal_var=False)
    lqrts.append(x[1])
    i+=1
    
    
normaliseddf['lqrt']=lqrts
'''   


normaliseddf['Wald p-value (Deseq2 Test)']=waldvalues

normaliseddf['mann-whitney u']=mannwhitneys



#correcting l2fc for volcano
normaliseddf['l2fc2']=normaliseddf['l2fc']-1




normaliseddf=normaliseddf.sort_values(by=['welch correction'])
normaliseddf

normaliseddf.to_csv('Simon_customdiffexp_final.csv')