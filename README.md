# BRWCP
BRCWP: Drug-Protein Interaction Prediction by Correcting the Effect of Incomplete Information in Heterogeneous Information
# Quick Start
to reproduce our results:
1. Download and unzip our code files
2. Run dti_devide.py to reproduce DTI data and their indexes for each fold and save them to ./dataset
3. Run run_main.py to reproduce the ten fold mean true positive rate list, false positive rate list, precision list, and reacll list about drug of the prediction results, and save them to ./result  
4. Run draw_fig.py to reproduce the AUC and AUPR of the prediction result and get the ROC/PRC firgure as follow
![Figure_1](https://user-images.githubusercontent.com/103353319/191415081-b9a80df1-c622-4dcf-9178-b1ca35a9ff84.png)
## parameters in run_main.py are:
\# rs: restart probability for random walk, defult=0.6 \
\# fr: dimensionals of low-dimensional feature of drugs, defult=500 \
\# fp: dimensionals of low-dimensional feature of proteins, defult=300 \
\# K: number of neighbors retained after pruning, defult=20 \
\# $\eta$: decay term of neighbors' weight, defult=0.7 \
\# $l_1$: number of random walks in the incomplete information network, defult=5\
\# $l_2$: number of random walks in the complete information network, defult=5

# Contacts
If you have any questions, please email Li Yanfei (lyfinf@163.com)
