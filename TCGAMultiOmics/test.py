from TCGAMultiOmics.multiomics import MultiOmicsData

folder_path = "/home/jonny_admin/PycharmProjects/nuclei-segmentation/data/tcga-assembler/LUAD"
external_data_path = "/home/jonny_admin/PycharmProjects/nuclei-segmentation/data/external"
luad_data = MultiOmicsData(cancer_type="LUAD", tcga_data_path=folder_path, external_data_path=external_data_path,
                           modalities=["GE", "MIR", "LNC", "CNV", "SNP", "PRO"])
