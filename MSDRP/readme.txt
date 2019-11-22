MSDRP: Modular within and between Score for Drug Response Prediction

Developer: Shiming Wang, Jie Li
           School of Computer Science and Technology, Harbin Institute of Technology, No. 92 Xidazhi Street, Nangang District, Harbin City, Heilongjiang Province, China, 150001. Email: shimingwang@hit.edu.cn, jieli@hit.edu.cn

This method is used to predict drug responses in cancer cell lines, based on a constructed heterogeneous drug-cell line network with multiple information. The network contains 251 drugs and 990 cell lines and 23774 associations.

This code is written by C++ and can be run on CodeBlocks 16.01 or later.

Input data:
1) "cell line-ID.txt": ID of 990 cell lines. 
2) "cell line-ID-962.txt": ID of 962 cell lines. These cell lines have gene expression profile.
3) "similarity_cellline-962.txt": The similarity of the 962 cell lines.
4) "drug_name.txt": 251 drug names.
5) "drug_name-187.txt": 187 drug names. These drugs have drug chemical structure feature.
6) "drug_name-with-target.txt": Name of drugs in the drug-target interaction network.
7) "drug-similarity-target-network.txt":  Drug similarity calculated based on the drug-target interaction network.
8) "similarity_drug-187.txt": Drug similarity calculated based on drug chemical structure feature.
9) "drug-cellline_association.txt":  The drug-cell line association network, including 251 drugs, 990 cell lines and 23774 associations.

Output data:
1) "unloocv.txt": The association score matrix of the drug-cell line association network.
2) "loocv.txt": The LOOCV result.

predict-association.xlsx: The predicted drug-cell line interactions. After LOOCV, choose the threshold when the Youden index reaches the maximum as the classification threshold. For the drug-cell line pair which is ¡°0¡± in the adjacency matrix, if its association score is larger than the threshold, it is a predicted drug-cell line interaction.
