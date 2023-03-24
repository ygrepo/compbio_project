import pandas as pd

gse_file = ''
gse = pd.read_csv(gse_file,index_col=0)
gse.head()

# strip symbol with version numbers.
index = []
for gene in gse.index.values:
    gene = str(gene)
    gene = gene.split('.')[0]
    index.append(gene)
gse.index = index

# collapse duplication
expression_dict = {}
running_length = len(gse)
for i in range(running_length):
    line = gse.iloc[i]
    gene = line.name
    try:
        gene_array = expression_dict[gene]
        gene_array = gene_array + line
        expression_dict[gene] = gene_array
    except:
        expression_dict[gene] = line
        
# transform dict to table
gene_list = []
array = []
for key,value in table.items():
    gene_list.append(key)
    array.append(value)
collapsed_table = pd.DataFrame(data=array,index=gene_list)
collapsed_table= collapsed_table.fillna(0)

collapsed_table.to_csv(gse_file)

