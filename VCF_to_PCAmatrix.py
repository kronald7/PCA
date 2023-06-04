### The final output of this py script is a matrix shaped csv file
# ready for the PCA analysis.

from pysam import VariantFile
import numpy as np
from sklearn import decomposition
import pandas as pd

vcf_filename = "/Users/ronky/Projects/WGS_SANRUHRP/snps_forPCA.vcf.gz"
panel_filename = "/Users/ronky/Projects/WGS_SANRUHRP/panel_file.txt" 

genotypes = []
samples = []
variant_ids = []
with VariantFile(vcf_filename) as vcf_reader:
    for record in vcf_reader:
        alleles = [record.samples[x].allele_indices for x in record.samples]
        samples = [sample for sample in record.samples]
        genotypes.append(alleles)
        variant_ids.append(record.id)

# make a label file with dictionary like {sample_id: population_code}

with open(panel_filename) as panel_file:
    labels = {}  # {sample_id: population_code}
    for line in panel_file:
        line = line.strip().split('\t')
        labels[line[0]] = line[1]

genotypes = np.array(genotypes)
matrix = np.count_nonzero(genotypes, axis=2)
matrix = matrix.T

pca = decomposition.PCA(n_components=2)
pca.fit(matrix)
print(pca.singular_values_)
to_plot = pca.transform(matrix)
print(to_plot.shape)

df = pd.DataFrame(matrix, columns=variant_ids, index=samples)
df['Population code'] = df.index.map(labels)
df.to_csv("pca_matrix.csv")