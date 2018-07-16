import os

import pandas as pd

TUMOR_NORMAL = 'tumor_normal'

CLINICAL_DRUG_FILENAME = "nationwidechildrens.org_clinical_drug.txt"
CLINICAL_PATIENT_FILENAME = "nationwidechildrens.org_clinical_patient.txt"

BCR_PATIENT_BARCODE = "bcr_patient_barcode"
HISTOLOGIC_SUBTYPE = "histologic_subtype"
PATHOLOGIC_STAGE = "pathologic_stage"


class ClinicalData:
    clinical_patient_colsname = ['bcr_patient_barcode', 'gender', 'race', 'histologic_diagnosis.1',
                                 'ajcc_pathologic_tumor_stage'
                                 ]

    pathologic_stage_map = {'Stage IA': 'Stage I', 'Stage IB': 'Stage I',
                            'Stage IIA': 'Stage II', 'Stage IIB': 'Stage II',
                            'Stage IIIA': 'Stage III', 'Stage IIIB': 'Stage III'}

    biospecimen_sample_colsname = ['bcr_sample_barcode', 'sample_type']

    clinical_drug_colsname = ['bcr_patient_barcode', 'pharmaceutical_therapy_drug_name', 'pharmaceutical_therapy_type',
                              'treatment_best_response']

    def __init__(self, cancer_type, folder_path):
        self.cancer_type = cancer_type

        # Import patients
        self.patient = pd.read_table(os.path.join(folder_path, "%s" % CLINICAL_PATIENT_FILENAME),
                                     sep="\t",
                                     skiprows=[1, 2],
                                     na_values=["[Not Available]", "[Not Applicable]"],
                                     usecols=ClinicalData.clinical_patient_colsname
                                     )
        self.patient.index = self.patient[BCR_PATIENT_BARCODE]
        self.patient.rename({"ajcc_pathologic_tumor_stage": ("%s" % PATHOLOGIC_STAGE),
                             "histologic_diagnosis.1": ("%s" % HISTOLOGIC_SUBTYPE)}, axis=1, inplace=True)
        self.patient.replace({('%s' % PATHOLOGIC_STAGE): ClinicalData.pathologic_stage_map}, inplace=True)


        # # Import biospecimen samples (not all samples included in dataset)
        # self.biospecimen = pd.read_table(os.path.join(folder_path, "genome.wustl.edu_biospecimen_sample.txt"),
        #                                  sep="\t",
        #                                  skiprows=[1, ],
        #                                  na_values="[Not Available]",
        #                                  usecols=ClinicalData.biospecimen_sample_colsname
        #                                  )
        # self.biospecimen.index = self.biospecimen["bcr_sample_barcode"]

        # Import clinical drug
        self.drugs = pd.read_table(os.path.join(folder_path, "%s" % CLINICAL_DRUG_FILENAME),
                                   sep="\t",
                                   skiprows=[1, 2],
                                   na_values=["[Not Available]", "[Unknown]", "[Not Applicable]"],
                                   usecols=ClinicalData.clinical_drug_colsname
                                   )
        self.drugs.index = self.drugs[BCR_PATIENT_BARCODE]

        # Save index's
        self.patient_barcodes = self.patient[BCR_PATIENT_BARCODE].tolist()
        # self.sample_barcodes = self.biospecimen["bcr_sample_barcode"].tolist()
        # self.drug_barcodes = self.biospecimen["bcr_sample_barcode"].tolist()

    def build_clinical_samples(self, all_samples):
        # Build table with samples clinical data from patients
        self.samples = pd.DataFrame(index=all_samples)

        self.samples[BCR_PATIENT_BARCODE] = self.samples.index.str[:-4]

        no_samples = self.samples.shape[0]

        # Merge patients clinical data with patient barcode as index
        # target = pd.merge(target, self.patient,
        #                      how="left", left_on="patient_barcode", right_on="patient_barcode")
        self.samples = self.samples.join(self.patient, on=BCR_PATIENT_BARCODE, how="left", rsuffix="_")
        if self.samples.shape[0] != no_samples:
            raise Exception("Clinical data merging has wrong number of samples")
        self.samples.drop(BCR_PATIENT_BARCODE+"_", axis=1, inplace=True)  # Remove redundant column
        # self.samples.dropna(axis=0, subset=["bcr_patient_barcode"], inplace=True) # Remove samples without clinical data

        self.samples = self.samples[self.samples[PATHOLOGIC_STAGE] != "[Discrepancy]"]
        self.samples.loc[self.samples.index.str.contains("-11"),
                         ('%s' % TUMOR_NORMAL)] = "Normal"  # Change stage label of normal samples to "Normal"
        self.samples.loc[self.samples.index.str.contains("-01"),
                         TUMOR_NORMAL] = "Tumor"  # Change stage label of normal samples to "Normal"


    def get_patient_barcodes(self):
        return self.patient_barcodes

    # def get_sample_barcodes(self):
    #     return self.sample_barcodes


