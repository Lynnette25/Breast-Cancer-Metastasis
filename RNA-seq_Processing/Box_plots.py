################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
# 	  1) Reads gene expression data from a CSV file.
# 	  2) Generates box plots for specified upregulated genes across different tumor types.
# 	  3) Performs statistical tests to compare expression levels between tumor types and adds significance bars to the plots.
# 	  4) Saves the generated plots to specified file paths.
################################################################################


#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# In[ ]:





# In[2]:


df=pd.read_csv("~/Metastasis/RNA/PE_counts/FULL_gene_data_PE_MET_GTEX_one_to_one.csv", index_col="gene_name")
df = df[~df.index.duplicated(keep='first')]
df.head()


# In[5]:


#pd.set_option("display.max_columns", None)
print(df.loc["C3orf80"]["Geneid"])


# In[11]:


import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

Up = ["KLF5", "KLF12", "SP3", "KLF10", "SP5","SP1", "SP4", "KLF14", "KLF15"]
colors = {
    'Liver': '#c3e8f8',   # Azul claro
    'Lung': '#cbf2d1',    # Verde Claro
    'Breast': '#ffbaca'   # Rosita
}

for TF in Up:
    # Extract tumor values
    liver_tumor_cols = [col for col in df.columns if 'Liver_Tumor' in col]
    Liver_Tumor_values = df.loc[TF, liver_tumor_cols].astype(float).tolist()

    lung_tumor_cols = [col for col in df.columns if 'Lung_Tumor' in col]
    Lung_Tumor_values = df.loc[TF, lung_tumor_cols].astype(float).tolist()

    breast_tumor_cols = [col for col in df.columns if 'Breast_Tumor' in col]
    Breast_Tumor_values = df.loc[TF, breast_tumor_cols].astype(float).tolist()

    # Optional: convert to int if all values are whole numbers
    Liver_Tumor_values = [int(val) if val.is_integer() else val for val in Liver_Tumor_values]
    Lung_Tumor_values = [int(val) if val.is_integer() else val for val in Lung_Tumor_values]
    Breast_Tumor_values = [int(val) if val.is_integer() else val for val in Breast_Tumor_values]

    data = [Liver_Tumor_values, Lung_Tumor_values, Breast_Tumor_values]
    labels = ['Liver', 'Lung', 'Breast']

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_axes([0, 0, 1, 1])

    # Create boxplot with custom edge colors
    bp = ax.boxplot(data, patch_artist=True)

    # Set custom border (edge) colors for each box
    for patch, color in zip(bp['boxes'], [colors['Liver'], colors['Lung'], colors['Breast']]):
        patch.set_facecolor('none')
        patch.set_edgecolor(color)
        patch.set_linewidth(2)


    ax.set_xticklabels(['Liver', 'Lung', 'Breast'])
    plt.title(f'Expression of {TF} in Tumor Types')
    plt.ylabel('Counts')

    # Optionally, add statistical significance
    p1 = stats.ttest_ind(Liver_Tumor_values, Lung_Tumor_values).pvalue
    p2 = stats.ttest_ind(Liver_Tumor_values, Breast_Tumor_values).pvalue
    p3 = stats.ttest_ind(Lung_Tumor_values, Breast_Tumor_values).pvalue

    def add_significance_bar(data1, data2, x1, x2, p_value):
        y_max = max(max(data1), max(data2))
        y = y_max * 1.1
        ax.plot([x1, x1, x2, x2], [y, y*1.02, y*1.02, y], color='black')
        ax.text((x1 + x2) / 2, y*1.04, f'p = {p_value:.4f}', ha='center', fontsize=8)

    add_significance_bar(Liver_Tumor_values, Lung_Tumor_values, 1, 2, p1)
    add_significance_bar(Liver_Tumor_values, Breast_Tumor_values, 1, 3, p2)
    add_significance_bar(Lung_Tumor_values, Breast_Tumor_values, 2, 3, p3)

    plt.savefig(f"~/Metastasis/Back_to_peaks/ {TF} in Liver or Lung.png", bbox_inches='tight')
    plt.show()


# In[8]:


TF = "CES1"

liver_tumor_cols = [col for col in df.columns if 'Liver_Tumor' in col]
Liver_Tumor_values = df.loc[TF, liver_tumor_cols].astype(float).tolist()

lung_tumor_cols = [col for col in df.columns if 'Lung_Tumor' in col]
Lung_Tumor_values = df.loc[TF, lung_tumor_cols].astype(float).tolist()

breast_tumor_cols = [col for col in df.columns if 'Breast_Tumor' in col]
Breast_Tumor_values = df.loc[TF, breast_tumor_cols].astype(float).tolist()

Normal_liver_cols = [col for col in df.columns if 'Normal_Liver' in col]
Normal_Liver_values = df.loc[TF, Normal_liver_cols].astype(float).tolist()

Normal_lung_cols = [col for col in df.columns if 'Normal_Lung' in col]
Normal_Lung_values = df.loc[TF, Normal_lung_cols].astype(float).tolist()

Normal_breast_cols = [col for col in df.columns if 'Normal_Breast' in col]
Normal_Breast_values = df.loc[TF, Normal_breast_cols].astype(float).tolist()

# Convert all lists to integers if they are whole numbers
Liver_Tumor_values = [int(val) if val.is_integer() else val for val in Liver_Tumor_values]
Lung_Tumor_values = [int(val) if val.is_integer() else val for val in Lung_Tumor_values]
Breast_Tumor_values = [int(val) if val.is_integer() else val for val in Breast_Tumor_values]
Normal_Liver_values = [int(val) if val.is_integer() else val for val in Normal_Liver_values]
Normal_Lung_values = [int(val) if val.is_integer() else val for val in Normal_Lung_values]
Normal_Breast_values = [int(val) if val.is_integer() else val for val in Normal_Breast_values]

# Print the values
print(Liver_Tumor_values)
print(Lung_Tumor_values)
print(Breast_Tumor_values)
print(Normal_Liver_values)
print(Normal_Lung_values)
print(Normal_Breast_values)



# In[5]:


# Extracting data for different groups


data = [Normal_Breast_values, Breast_Tumor_values, Normal_Liver_values, Liver_Tumor_values,Normal_Lung_values, Lung_Tumor_values]

fig = plt.figure(figsize=(6, 2.4))

# Creating axes instance
ax = fig.add_axes([0, 0, 1, 1])

# Creating box plot
bp = ax.boxplot(data)
ax.set_xticklabels([ 'Normal Breast', 'Breast Tumor', 'Normal Liver', 'Liver Tumor' , 'Normal Lung', 'Lung Tumor'])
#ax.set_ylim(800, 1020)


plt.title(f'Box Plot of {TF} Expression in Metastatic Cancer')
plt.ylabel('Counts')

# Add significance bars
def add_significance_bar(data1, data2, v, x1, x2, p_value):
    """Add significance bar between two boxplot groups."""
    # Define the height of the significance bar
    data_max = max(max(data1), max(data2))
    y1 = data_max + v
    y2 = data_max + v
    h = 0.06 * (y1 - min(data1 + data2))  # Adjust height based on overall data range
    #h =100

    # Draw the bar
    ax.plot([x1, x1, x2, x2], [y1 - h, y1, y2, y2 - h], color='black')
    # Annotate the p-value
    ax.text((x1 + x2) / 2, y2+6 , f'p = {p_value:.4f}', ha='center')

# Perform statistical tests and add significance bars
p_value_1 = stats.ttest_ind(Normal_Breast_values, Breast_Tumor_values).pvalue
p_value_2 = stats.ttest_ind( Normal_Liver_values, Liver_Tumor_values).pvalue
p_value_3 = stats.ttest_ind(Normal_Lung_values, Lung_Tumor_values).pvalue
add_significance_bar(Normal_Breast_values, Breast_Tumor_values, 150, 1, 2, p_value_1)
add_significance_bar(Normal_Liver_values, Liver_Tumor_values, 150, 3, 4, p_value_2)
add_significance_bar(Normal_Lung_values, Lung_Tumor_values, 150, 5,6, p_value_3)

plt.savefig(f"/Users/francesheredia/Documents/Metastasis/KM_Down_Filtered/Box Plot of {TF} A in Metastatic Cancer.png", bbox_inches='tight')
plt.show()


# In[114]:


Up = ["C3orf80"]
for TF in Up:

    liver_tumor_cols = [col for col in df.columns if 'Liver_Tumor' in col]
    Liver_Tumor_values = df.loc[TF, liver_tumor_cols].astype(float).tolist()

    lung_tumor_cols = [col for col in df.columns if 'Lung_Tumor' in col]
    Lung_Tumor_values = df.loc[TF, lung_tumor_cols].astype(float).tolist()

    breast_tumor_cols = [col for col in df.columns if 'Breast_Tumor' in col]
    Breast_Tumor_values = df.loc[TF, breast_tumor_cols].astype(float).tolist()

    Normal_liver_cols = [col for col in df.columns if 'Normal_Liver' in col]
    Normal_Liver_values = df.loc[TF, Normal_liver_cols].astype(float).tolist()

    Normal_lung_cols = [col for col in df.columns if 'Normal_Lung' in col]
    Normal_Lung_values = df.loc[TF, Normal_lung_cols].astype(float).tolist()

    Normal_breast_cols = [col for col in df.columns if 'Normal_Breast' in col]
    Normal_Breast_values = df.loc[TF, Normal_breast_cols].astype(float).tolist()

    # Convert all lists to integers if they are whole numbers
    Liver_Tumor_values = [int(val) if val.is_integer() else val for val in Liver_Tumor_values]
    Lung_Tumor_values = [int(val) if val.is_integer() else val for val in Lung_Tumor_values]
    Breast_Tumor_values = [int(val) if val.is_integer() else val for val in Breast_Tumor_values]
    Normal_Liver_values = [int(val) if val.is_integer() else val for val in Normal_Liver_values]
    Normal_Lung_values = [int(val) if val.is_integer() else val for val in Normal_Lung_values]
    Normal_Breast_values = [int(val) if val.is_integer() else val for val in Normal_Breast_values]

    data = [Normal_Breast_values, Breast_Tumor_values, Normal_Liver_values, Liver_Tumor_values,Normal_Lung_values, Lung_Tumor_values]

    fig = plt.figure(figsize=(6, 2.4))

    # Creating axes instance
    ax = fig.add_axes([0, 0, 1, 1])

    # Creating box plot
    bp = ax.boxplot(data)
    ax.set_xticklabels([ 'Normal Breast', 'Breast Tumor', 'Normal Liver', 'Liver Tumor' , 'Normal Lung', 'Lung Tumor'])
    #ax.set_ylim(800, 1020)


    plt.title(f'Box Plot of {TF} Expression in Metastatic Cancer')
    plt.ylabel('Counts')

    # Add significance bars
    def add_significance_bar(data1, data2, v, x1, x2, p_value):
        """Add significance bar between two boxplot groups."""
        # Define the height of the significance bar
        data_max = max(max(data1), max(data2))
        y1 = data_max * 1.20 + v
        y2 = data_max * 1.20 + v
        h = 0.06 * (y1 - min(data1 + data2))  # Adjust height based on overall data range
        #h =100

        # Draw the bar
        ax.plot([x1, x1, x2, x2], [y1 - h, y1*1.1, y2*1.1, y2 - h], color='black')
        # Annotate the p-value
        ax.text((x1 + x2) / 2, data_max * 1.20 *1.1 , f'p = {p_value:.4f}', ha='center')

    # Perform statistical tests and add significance bars
    p_value_1 = stats.ttest_ind(Normal_Breast_values, Breast_Tumor_values).pvalue
    p_value_2 = stats.ttest_ind( Normal_Liver_values, Liver_Tumor_values).pvalue
    p_value_3 = stats.ttest_ind(Normal_Lung_values, Lung_Tumor_values).pvalue
    add_significance_bar(Normal_Breast_values, Breast_Tumor_values, 0, 1, 2, p_value_1)
    add_significance_bar(Normal_Liver_values, Liver_Tumor_values, 0, 3, 4, p_value_2)
    add_significance_bar(Normal_Lung_values, Lung_Tumor_values, 0, 5,6, p_value_3)

    plt.savefig(f"~/Metastasis/RNA/PE_counts/Box_Plots/Box Plot of {TF}  in Metastic Liver Specific Genes.png", bbox_inches='tight')
    plt.show()




# In[6]:


liver = pd.read_csv("/Users/francolab3/Documents/Metastasis/RNA/PE_counts/KM_Liver_UP_Filtered_cox_results_PE.csv")
liver.head()


# In[7]:


Up = liver["TF"]
for TF in Up:

    liver_tumor_cols = [col for col in df.columns if 'Liver_Tumor' in col]


    Liver_Tumor_values  = df.loc[TF, liver_tumor_cols].tolist()

    lung_tumor_cols = [col for col in df.columns if 'Lung_Tumor' in col]


    Lung_Tumor_values  = df.loc[TF, lung_tumor_cols].tolist()

    Breast_tumor_cols = [col for col in df.columns if 'Breast_Tumor' in col]


    Breast_Tumor_values  = df.loc[TF, Breast_tumor_cols].tolist()



    Normal_liver_cols = [col for col in df.columns if 'Normal_Liver' in col]


    Normal_Liver_values  = df.loc[TF, Normal_liver_cols].tolist()

    Normal_lung_cols = [col for col in df.columns if 'Normal_Lung' in col]


    Normal_Lung_values  = df.loc[TF, Normal_lung_cols].tolist()

    Normal_breast_cols = [col for col in df.columns if 'Normal_Breast' in col]


    Normal_Breast_values  = df.loc[TF, Normal_breast_cols ].tolist()


    data = [Normal_Breast_values, Breast_Tumor_values, Normal_Liver_values, Liver_Tumor_values,Normal_Lung_values, Lung_Tumor_values]

    fig = plt.figure(figsize=(6, 2.4))

    # Creating axes instance
    ax = fig.add_axes([0, 0, 1, 1])

    # Creating box plot
    bp = ax.boxplot(data)
    ax.set_xticklabels([ 'Normal Breast', 'Breast Tumor', 'Normal Liver', 'Liver Tumor' , 'Normal Lung', 'Lung Tumor'])
    #ax.set_ylim(800, 1020)


    plt.title(f'Box Plot of {TF} Expression in Metastatic Cancer')
    plt.ylabel('Counts')

    # Add significance bars
    def add_significance_bar(data1, data2, v, x1, x2, p_value):
        """Add significance bar between two boxplot groups."""
        # Define the height of the significance bar
        data_max = max(max(data1), max(data2))
        y1 = data_max * 1.20 + v
        y2 = data_max * 1.20 + v
        h = 0.06 * (y1 - min(data1 + data2))  # Adjust height based on overall data range
        #h =100

        # Draw the bar
        ax.plot([x1, x1, x2, x2], [y1 - h, y1*1.1, y2*1.1, y2 - h], color='black')
        # Annotate the p-value
        ax.text((x1 + x2) / 2, data_max * 1.20 *1.1 , f'p = {p_value:.4f}', ha='center')

    # Perform statistical tests and add significance bars
    p_value_1 = stats.ttest_ind(Normal_Breast_values, Breast_Tumor_values).pvalue
    p_value_2 = stats.ttest_ind( Normal_Liver_values, Liver_Tumor_values).pvalue
    p_value_3 = stats.ttest_ind(Normal_Lung_values, Lung_Tumor_values).pvalue
    add_significance_bar(Normal_Breast_values, Breast_Tumor_values, 0, 1, 2, p_value_1)
    add_significance_bar(Normal_Liver_values, Liver_Tumor_values, 0, 3, 4, p_value_2)
    add_significance_bar(Normal_Lung_values, Lung_Tumor_values, 0, 5,6, p_value_3)

    plt.savefig(f"/Users/francolab3/Documents/Metastasis/RNA/PE_counts/Box_Plots/Box Plot of {TF}  in Metastic Liver Specific Genes.png", bbox_inches='tight')
    plt.show()




# In[11]:


MET = pd.read_csv("~/Metastasis/RNA/PE_counts/KM_MET_UP_Filtered_cox_results_PE.csv")

Up = MET["TF"]
for TF in Up:

    liver_tumor_cols = [col for col in df.columns if 'Liver_Tumor' in col]


    Liver_Tumor_values  = df.loc[TF, liver_tumor_cols].tolist()

    lung_tumor_cols = [col for col in df.columns if 'Lung_Tumor' in col]


    Lung_Tumor_values  = df.loc[TF, lung_tumor_cols].tolist()

    Breast_tumor_cols = [col for col in df.columns if 'Breast_Tumor' in col]


    Breast_Tumor_values  = df.loc[TF, Breast_tumor_cols].tolist()



    Normal_liver_cols = [col for col in df.columns if 'Normal_Liver' in col]


    Normal_Liver_values  = df.loc[TF, Normal_liver_cols].tolist()

    Normal_lung_cols = [col for col in df.columns if 'Normal_Lung' in col]


    Normal_Lung_values  = df.loc[TF, Normal_lung_cols].tolist()

    Normal_breast_cols = [col for col in df.columns if 'Normal_Breast' in col]


    Normal_Breast_values  = df.loc[TF, Normal_breast_cols ].tolist()


    data = [Normal_Breast_values, Breast_Tumor_values, Normal_Liver_values, Liver_Tumor_values,Normal_Lung_values, Lung_Tumor_values]

    fig = plt.figure(figsize=(6, 2.4))

    # Creating axes instance
    ax = fig.add_axes([0, 0, 1, 1])

    # Creating box plot
    bp = ax.boxplot(data)
    ax.set_xticklabels([ 'Normal Breast', 'Breast Tumor', 'Normal Liver', 'Liver Tumor' , 'Normal Lung', 'Lung Tumor'])
    #ax.set_ylim(800, 1020)


    plt.title(f'Box Plot of {TF} Expression in Metastatic Cancer')
    plt.ylabel('Counts')

    # Add significance bars
    def add_significance_bar(data1, data2, v, x1, x2, p_value):
        """Add significance bar between two boxplot groups."""
        # Define the height of the significance bar
        data_max = max(max(data1), max(data2))
        y1 = data_max * 1.20 + v
        y2 = data_max * 1.20 + v
        h = 0.06 * (y1 - min(data1 + data2))  # Adjust height based on overall data range
        #h =100

        # Draw the bar
        ax.plot([x1, x1, x2, x2], [y1 - h, y1*1.1, y2*1.1, y2 - h], color='black')
        # Annotate the p-value
        ax.text((x1 + x2) / 2, data_max * 1.20 *1.1 , f'p = {p_value:.4f}', ha='center')

    # Perform statistical tests and add significance bars
    p_value_1 = stats.ttest_ind(Normal_Breast_values, Breast_Tumor_values).pvalue
    p_value_2 = stats.ttest_ind( Normal_Liver_values, Liver_Tumor_values).pvalue
    p_value_3 = stats.ttest_ind(Normal_Lung_values, Lung_Tumor_values).pvalue
    add_significance_bar(Normal_Breast_values, Breast_Tumor_values, 0, 1, 2, p_value_1)
    add_significance_bar(Normal_Liver_values, Liver_Tumor_values, 0, 3, 4, p_value_2)
    add_significance_bar(Normal_Lung_values, Lung_Tumor_values, 0, 5,6, p_value_3)

    plt.savefig(f"/Users/francolab3/Documents/Metastasis/RNA/PE_counts/Box_Plots/Box Plot of {TF}  in Metastic Liver Specific Genes.png", bbox_inches='tight')
    plt.show()




# In[111]:


lung = pd.read_csv("/Users/francolab3/Documents/Metastasis/RNA/PE_counts/ER_positive_Metabric_Cox_results.csv")
Up = lung["TF"]
for TF in Up:

    liver_tumor_cols = [col for col in df.columns if 'Liver_Tumor' in col]


    Liver_Tumor_values  = df.loc[TF, liver_tumor_cols].tolist()

    lung_tumor_cols = [col for col in df.columns if 'Lung_Tumor' in col]


    Lung_Tumor_values  = df.loc[TF, lung_tumor_cols].tolist()

    Breast_tumor_cols = [col for col in df.columns if 'Breast_Tumor' in col]


    Breast_Tumor_values  = df.loc[TF, Breast_tumor_cols].tolist()



    Normal_liver_cols = [col for col in df.columns if 'Normal_Liver' in col]


    Normal_Liver_values  = df.loc[TF, Normal_liver_cols].tolist()

    Normal_lung_cols = [col for col in df.columns if 'Normal_Lung' in col]


    Normal_Lung_values  = df.loc[TF, Normal_lung_cols].tolist()

    Normal_breast_cols = [col for col in df.columns if 'Normal_Breast' in col]


    Normal_Breast_values  = df.loc[TF, Normal_breast_cols ].tolist()


    data = [Normal_Breast_values, Breast_Tumor_values, Normal_Liver_values, Liver_Tumor_values,Normal_Lung_values, Lung_Tumor_values]

    fig = plt.figure(figsize=(6, 2.4))

    # Creating axes instance
    ax = fig.add_axes([0, 0, 1, 1])

    # Creating box plot
    bp = ax.boxplot(data)
    ax.set_xticklabels([ 'Normal Breast', 'Breast Tumor', 'Normal Liver', 'Liver Tumor' , 'Normal Lung', 'Lung Tumor'])
    #ax.set_ylim(800, 1020)


    plt.title(f'Box Plot of {TF} Expression in Metastatic Cancer')
    plt.ylabel('Counts')

    # Add significance bars
    def add_significance_bar(data1, data2, v, x1, x2, p_value):
        """Add significance bar between two boxplot groups."""
        # Define the height of the significance bar
        data_max = max(max(data1), max(data2))
        y1 = data_max * 1.20 + v
        y2 = data_max * 1.20 + v
        h = 0.06 * (y1 - min(data1 + data2))  # Adjust height based on overall data range
        #h =100

        # Draw the bar
        ax.plot([x1, x1, x2, x2], [y1 - h, y1*1.1, y2*1.1, y2 - h], color='black')
        # Annotate the p-value
        ax.text((x1 + x2) / 2, data_max * 1.20 *1.1 , f'p = {p_value:.4f}', ha='center')

    # Perform statistical tests and add significance bars
    p_value_1 = stats.ttest_ind(Normal_Breast_values, Breast_Tumor_values).pvalue
    p_value_2 = stats.ttest_ind( Normal_Liver_values, Liver_Tumor_values).pvalue
    p_value_3 = stats.ttest_ind(Normal_Lung_values, Lung_Tumor_values).pvalue
    add_significance_bar(Normal_Breast_values, Breast_Tumor_values, 0, 1, 2, p_value_1)
    add_significance_bar(Normal_Liver_values, Liver_Tumor_values, 0, 3, 4, p_value_2)
    add_significance_bar(Normal_Lung_values, Lung_Tumor_values, 0, 5,6, p_value_3)

    plt.savefig(f"/Users/francolab3/Documents/Metastasis/RNA/PE_counts/Box_Plots/Box Plot of {TF}  in Metastic Liver Specific Genes.png", bbox_inches='tight')
    plt.show()




# In[1]:


'''
Down =['MYD88' , 'ITGAX' , 'MALT1' , 'PARP14' , 'EVI2B' , 'RASGRP4' , 'EIF1B' , 'DENND5A' , 'KLHL6' , 'NCF2' , 'COL4A3' , 'VNN2' , 'GPCPD1' , 'DNAJB5' , 'FGR' , 'RAB8B' , 'FERMT3' , 'RASSF2' , 'PIK3R5' , 'SASH3' , 'MBNL1' , 'CD53' , 'BCL2A1' , 'CASP4' , 'BTN2A2' , 'CD300C' , 'RFTN1' , 'SH2B3' , 'INPP5D' , 'BIN2' , 'CXCL16' , 'LCP1' , 'STON2' , 'DOCK2' , 'PRKCB' , 'NAPSB' , 'FMNL1' , 'RBMS2' , 'HCLS1' , 'ZNF267' , 'LCP2' , 'MYO1F' , 'CRTAM' , 'ALOX5' , 'FCER1G' , 'TRIM21' , 'CTSS' , 'CARD16' , 'WDFY4' , 'AKR1B1' , 'SPI1' , 'LSAMP' , 'IL12RB1' , 'RGS18' , 'SEMA3B' , 'KLRD1' , 'NIPAL4' , 'CHST11' , 'CD86' , 'HCK' , 'BTK' , 'NFKBIA' , 'TPST2' , 'PIM2' , 'SIGLEC14' , 'ARHGAP26' , 'DOK3' , 'ACSL4' , 'TYROBP' , 'PLEKHO1' , 'CD101' , 'SNX20' , 'CXCR6' , 'SYT11' , 'NFE2' , 'LIF' , 'PTGIR' , 'SYNPO2' , 'PALLD' , 'ALOX5AP' , 'FYN' , 'CYTIP' , 'VENTX' , 'LST1' , 'ADAMTS8' , 'CCNJL' , 'MIAT' , 'RAI14' , 'CD274' , 'PTPN22' , 'CLIC5' , 'RAB33A' , 'SPHK1' , 'CLEC5A' , 'GFI1' , 'GFPT2' , 'PRDM1' , 'FPR1' , 'RELB' , 'NCF4' , 'CXCL5' , 'TGFB2' , 'JUNB' , 'ZEB1' , 'SLA' , 'KLHL29' , 'HYAL2' , 'SRGN' , 'ITGAM' , 'SLC16A5' , 'TNFAIP6' , 'ARHGAP9' , 'LATS2' , 'TACC1' , 'FCGR1B' , 'FES' , 'ADPRH' , 'ETV1' , 'PDE4D' , 'CD8A' , 'LY86' , 'BIRC3' , 'GBP2' , 'TBC1D10C' , 'MBIP' , 'RNASE2' , 'GBP5' , 'ITPRIP' , 'SAMSN1' , 'SH2D1B' , 'C1orf116' , 'CLEC12A' , 'GZMM' , 'MSC' , 'NFAM1' , 'SLC2A3' , 'ITK' , 'LRRK2' , 'CCL5' , 'NR4A3' , 'ARID5A' , 'COL17A1' , 'EDN1' , 'ZBP1' , 'ATP1A2' , 'NR4A1' , 'GADD45A' , 'SGCD' , 'CRISPLD2' , 'KCNMB4' , 'CAMK2A' , 'CCR4' , 'CTLA4' , 'TOX' , 'PLAT' , 'CLDN18' , 'PIK3CD' , 'PRPH2' , 'GJA5' , 'FGF7' , 'GPA33' , 'NEFM' , 'IRF1' , 'CD5' , 'LMO3' , 'MMP2' , 'PRELP' , 'ADAM6' , 'KLRC2' , 'SLC7A4' , 'TTC24' , 'ADAMTS4' , 'CCL7' , 'GPR31']

for TF in Down:

    liver_tumor_cols = [col for col in df.columns if 'Liver_Tumor' in col]


    Liver_Tumor_values  = df.loc[TF, liver_tumor_cols].tolist()

    lung_tumor_cols = [col for col in df.columns if 'Lung_Tumor' in col]


    Lung_Tumor_values  = df.loc[TF, lung_tumor_cols].tolist()

    Breast_tumor_cols = [col for col in df.columns if 'Breast_Tumor' in col]


    Breast_Tumor_values  = df.loc[TF, Breast_tumor_cols].tolist()



    Normal_liver_cols = [col for col in df.columns if 'Normal_Liver' in col]


    Normal_Liver_values  = df.loc[TF, Normal_liver_cols].tolist()

    Normal_lung_cols = [col for col in df.columns if 'Normal_Lung' in col]


    Normal_Lung_values  = df.loc[TF, Normal_lung_cols].tolist()

    Normal_breast_cols = [col for col in df.columns if 'Normal_Breast' in col]


    Normal_Breast_values  = df.loc[TF, Normal_breast_cols ].tolist()


    data = [Normal_Breast_values, Breast_Tumor_values, Normal_Liver_values, Liver_Tumor_values,Normal_Lung_values, Lung_Tumor_values]

    fig = plt.figure(figsize=(6, 2.4))

    # Creating axes instance
    ax = fig.add_axes([0, 0, 1, 1])

    # Creating box plot
    bp = ax.boxplot(data)
    ax.set_xticklabels([ 'Normal Breast', 'Breast Tumor', 'Normal Liver', 'Liver Tumor' , 'Normal Lung', 'Lung Tumor'])
    #ax.set_ylim(800, 1020)


    plt.title(f'Box Plot of {TF} Expression in Metastatic Cancer')
    plt.ylabel('Counts')

    # Add significance bars
    def add_significance_bar(data1, data2, v, x1, x2, p_value):
        """Add significance bar between two boxplot groups."""
        # Define the height of the significance bar
        data_max = max(max(data1), max(data2))
        y1 = data_max * 1.20 + v
        y2 = data_max * 1.20 + v
        h = 0.06 * (y1 - min(data1 + data2))  # Adjust height based on overall data range
        #h =100

        # Draw the bar
        ax.plot([x1, x1, x2, x2], [y1 - h, y1*1.1, y2*1.1, y2 - h], color='black')
        # Annotate the p-value
        ax.text((x1 + x2) / 2, data_max * 1.20 *1.1 , f'p = {p_value:.4f}', ha='center')

    # Perform statistical tests and add significance bars
    p_value_1 = stats.ttest_ind(Normal_Breast_values, Breast_Tumor_values).pvalue
    p_value_2 = stats.ttest_ind( Normal_Liver_values, Liver_Tumor_values).pvalue
    p_value_3 = stats.ttest_ind(Normal_Lung_values, Lung_Tumor_values).pvalue
    add_significance_bar(Normal_Breast_values, Breast_Tumor_values, 0, 1, 2, p_value_1)
    add_significance_bar(Normal_Liver_values, Liver_Tumor_values, 0, 3, 4, p_value_2)
    add_significance_bar(Normal_Lung_values, Lung_Tumor_values, 0, 5,6, p_value_3)

    plt.savefig(f"/Users/francolab3/Documents/Metastasis/KM_Liver_Down_Filtered/Box Plot of {TF} in Metastic Liver Specific Genes.png", bbox_inches='tight')
    plt.show()

    '''


# In[ ]:




