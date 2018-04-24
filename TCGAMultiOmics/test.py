from TCGAMultiOmics.multiomics import MultiOmicsData

folder_path = "/home/jonny_admin/PycharmProjects/nuclei-segmentation/data/tcga-assembler/LUAD"
luad_data = MultiOmicsData(cancer_type="LUAD", folder_path=folder_path,
                           modalities=["GE", "MIR", "LNC", "CNV", "SNP", "PRO"])
