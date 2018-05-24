from TCGAMultiOmics.multiomics import MultiOmicsData
from TCGAMultiOmics.network import HeterogeneousNetwork

folder_path = "/home/jonny_admin/PycharmProjects/nuclei-segmentation/data/tcga-assembler/LUAD"
external_data_path = "/home/jonny_admin/PycharmProjects/nuclei-segmentation/data/external"
luad_data = MultiOmicsData(cancer_type="LUAD", tcga_data_path=folder_path, external_data_path=external_data_path,
                           modalities=["GE",
                                       "MIR",
                                       # "LNC",
                                       # "CNV",
                                       # "SNP",
                                       # "PRO"
                                       ])

network = HeterogeneousNetwork(modalities=["MIR", "GE"], multi_omics=luad_data)
network.add_edges_from_modality(modality="GE")
print(len(network.G.nodes()))
print(len(network.G.edges()))
# print(luad_data.LNC.get_genes_list())