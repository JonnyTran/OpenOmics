from TCGAMultiOmics.multiomics import MultiOmicsData

folder_path = "/home/jonny_admin/PycharmProjects/Bioinformatics_ExternalData/tcga-assembler/LUAD"
external_data_path = "/home/jonny_admin/PycharmProjects/Bioinformatics_ExternalData/"
luad_data = MultiOmicsData(cancer_type="LUAD", tcga_data_path=folder_path, external_data_path=external_data_path,
                           process_genes_info=False,
                           modalities=[
                                       # "GE",
                                       # "MIR",
                                       "LNC",
                                       # "CNV",
                                       # "SNP",
                                       # "PRO"
                                       ])
# print(luad_data.load_data(modalities=["GE", "MIR", "LNC"]))
# print(len(luad_data.LNC.get_genes_list()))
print(luad_data.LNC.get_genes_info().head())
# print(luad_data.MIR.get_genes_info())
# print(luad_data.GE.get_genes_info().T.apply(lambda x: x.nunique(), axis=1))
