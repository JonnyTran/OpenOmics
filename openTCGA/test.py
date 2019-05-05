from openTCGA.multiomic import MultiOmicsData

folder_path = "/data/datasets/Bioinformatics_ExternalData/tcga-assembler/LUAD"
external_data_path = "/data/datasets/Bioinformatics_ExternalData/"
luad_data = MultiOmicsData(cancer_type="LUAD", tcga_data_path=folder_path, external_data_path=external_data_path,
                           process_genes_info=True,
                           modalities=[
                                       # "GE",
                                       "MIR",
                                       # "LNC",
                                       # "CNV",
                                       # "SNP",
                                       # "PRO"
                                       ])

# LNC = luad_data.LNC.get_genes_info()
# print(luad_data.load_data(modalities=["GE", "MIR", "LNC"]))
# print(len(luad_data.LNC.get_genes_list()))
# print(LNC.columns)
# print("LNC transcripts matched", LNC["Rfams"].notnull().sum())
# print(luad_data.MIR.get_genes_info())
# print(luad_data.GE.get_genes_info().T.apply(lambda x: x.nunique(), axis=1))

# table = pd.read_table(luad_data.LNC.lncBase_interactions_file_path)
# print("matching geneName", len(set(LNC.index) & set(table["geneName"])))
# print("matching gene_id", len(set(LNC.index) & set(table["geneId"])))

# print(luad_data.LNC.get_lncRInter_interactions())

# print(LNC.head())