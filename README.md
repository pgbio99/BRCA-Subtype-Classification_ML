Machine Learning and Bioinformatics Framework for Breast Cancer Subtype Classification and Biomarker Discovery

📌 Overview

This project presents an integrative approach combining **bioinformatics analysis** and **machine learning models** to identify reliable biomarkers and classify molecular subtypes of **breast cancer** (BRCA) using **microarray gene expression data** from GEO datasets.
By leveraging tools such as differential expression analysis, protein-protein interaction networks, survival analysis, and supervised machine learning (Random Forest & KNN), the study uncovers 35 potential subtype-specific biomarkers and evaluates their diagnostic and prognostic significance.

🔬 Objectives

Identify DEGs linked to breast cancer subtypes using GEO microarray data, then construct PPI networks and identify key hub genes.
Perform subtype-specific DEG analysis to find unique molecular signatures, followed by developing ML models for subtype classification and biomarker prioritization.

🧪 Datasets Used

**GSE86374** (10 normal, 50 tumor samples)
**GSE57297** (7 normal, 25 tumor samples)  
 (Retrieved from the Gene Expression Omnibus (GEO))

🧰 Technologies and Tools

 **Languages**: Python, R  
 **Libraries**: `scikit-learn`, `matplotlib`, `pandas`, `seaborn`, `Limma` (R)  
 **Databases**: GEO, STRING, Enrichr, TNMplot, GEPIA2, cBioPortal 
 **Software**: Cytoscape, Google Colab

🧬 Bioinformatic Analysis

 **72 common DEGs** identified between datasets
 **10 hub genes** highlighted from the PPI network (e.g., KIT, TPM3, MYLK, COL10A1)
 **Survival analysis** confirmed clinical significance of multiple genes.
  
🧠 Machine Learning Models

| Model         | Accuracy | F1 Score | Sensitivity | Specificity |
|---------------|----------|----------|-------------|-------------|
| Random Forest | ~91%     | ~89%     | ~90%        | ~92%        |
| KNN           | ~86%     | ~84%     | ~85%        | ~87%        |

- **Random Forest** showed superior performance on microarray data.
- Feature importance was used to prioritize candidate biomarker genes.

**5276 subtype specific DEGs** involved in ML modeling
 **35 candidate biomarkers** prioritized via ML, including:  
  `TMEM45B`, `POF1B`, `STIL`, `DNA2`, `E2F7`, `RARA`, `HSPA5`, etc.

 📂 Project Structure

 BRCA-Subtype-Classification/
│
├── README.md                          
├── Random_Forest_BRCA_Subtypes.ipynb 
├── Project_Report_Prashansha_Goel.pdf
│
├── data/                              
    ├── GSE86374_normalized.csv
    └── GSE57297_normalized.csv
