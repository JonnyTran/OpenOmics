import os

import pandas as pd

TUMOR = "Tumor"
NORMAL = "Normal"

BCR_PATIENT_BARCODE = "bcr_patient_barcode"
HISTOLOGIC_SUBTYPE = "histologic_subtype"
PATHOLOGIC_STAGE = "pathologic_stage"
TUMOR_NORMAL = 'tumor_normal'
PREDICTED_SUBTYPE = 'predicted_subtype'

class ClinicalData:
    pathologic_stage_map = {'Stage IA': 'Stage I', 'Stage IB': 'Stage I',
                            'Stage IIA': 'Stage II', 'Stage IIB': 'Stage II',
                            'Stage IIIA': 'Stage III', 'Stage IIIB': 'Stage III'}

    def __init__(self, cohort_name, patients_file, patient_column="bcr_patient_barcode",
                 columns=None):
        """

        Args:
            cohort_name:
            folder_path:
            patients_file:
            patient_column:
            columns (list): default None. Example: ['bcr_patient_barcode', 'gender', 'race', 'histologic_diagnosis', 'tumor_status', 'death_days_to',
                          'ajcc_pathologic_tumor_stage']
        """
        self.cohort_name = cohort_name
        self.patient_column = patient_column

        if not os.path.exists(patients_file):
            raise FileNotFoundError(patients_file)

        self.patient = pd.read_table(patients_file,
                                     sep="\t",
                                     skiprows=[1, 2],
                                     na_values=["[Not Available]", "[Unknown]", "[Not Applicable]", "[Discrepancy]"],
                                     usecols=columns
                                     )

        self.patient_barcodes = self.patient[patient_column].tolist()
        self.patient.set_index(patient_column, inplace=True)

        # Rename columns
        self.patient.rename({"ajcc_pathologic_tumor_stage": ("%s" % PATHOLOGIC_STAGE),
                             "histological_type": ("%s" % HISTOLOGIC_SUBTYPE),
                             "histologic_diagnosis.1": ("%s" % HISTOLOGIC_SUBTYPE)}, axis=1, inplace=True)
        self.patient.replace({PATHOLOGIC_STAGE: ClinicalData.pathologic_stage_map}, inplace=True)

    @classmethod
    def name(self):
        return self.__class__.__name__

    def build_clinical_samples(self, all_samples, index="bcr_patient_barcode"):
        # Build table with samples clinical data from patients
        self.samples = pd.DataFrame(index=all_samples)

        self.samples[index] = self.samples.index.str[:-4]

        no_samples = self.samples.shape[0]

        # Merge patients clinical data with patient barcode as index
        # target = pd.merge(target, self.patient,
        #                      how="left", left_on="patient_barcode", right_on="patient_barcode")
        self.samples = self.samples.join(self.patient, on=index, how="left", rsuffix="_")
        if self.samples.shape[0] != no_samples:
            raise Exception("Clinical data merging has wrong number of samples")
        self.samples.drop(BCR_PATIENT_BARCODE+"_", axis=1, inplace=True)  # Remove redundant column
        # self.samples.dropna(axis=0, subset=["bcr_patient_barcode"], inplace=True) # Remove samples without clinical data

        self.samples = self.samples[self.samples[PATHOLOGIC_STAGE] != "[Discrepancy]"]
        self.samples.loc[self.samples.index.str.contains("-11"),
                         ('%s' % TUMOR_NORMAL)] = NORMAL  # Change stage label of normal samples to "Normal"
        self.samples.loc[self.samples.index.str.contains("-01"),
                         TUMOR_NORMAL] = TUMOR  # Change stage label of normal samples to "Normal"

    def add_drug_response_data(self, file_path="nationwidechildrens.org_clinical_drug.txt", patient_column="bcr_patient_barcode",
                               columns=['bcr_patient_barcode', 'pharmaceutical_therapy_drug_name', 'pharmaceutical_therapy_type', 'treatment_best_response'],
                               drug_name_col=None, response_column=None):
        if not os.path.exists(file_path):
            raise FileNotFoundError(file_path)

        self.drug_name_col = drug_name_col
        self.response_column = response_column

        self.drugs = pd.read_table(file_path,
                                   sep="\t",
                                   skiprows=[1, 2],
                                   na_values=["[Not Available]", "[Unknown]", "[Not Applicable]"],
                                   usecols=columns
                                   )
        self.drugs.set_index(patient_column, inplace=True)

    def add_biospecimen_data(self, file_path="genome.wustl.edu_biospecimen_sample.txt", patient_col_name="bcr_patient_barcode",
                             columns=['bcr_sample_barcode', 'sample_type']):
        if not os.path.exists(file_path):
            raise FileNotFoundError(file_path)

        self.biospecimen = pd.read_table(file_path, sep="\t", skiprows=[1, ],
                                         na_values=["[Not Available]", "[Unknown]", "[Not Applicable]"],
                                         usecols=columns
                                         )
        self.sample_barcodes = self.biospecimen[patient_col_name].tolist()
        self.biospecimen.set_index(patient_col_name, inplace=True)


    def get_patient_barcodes(self):
        return self.patient_barcodes

    def get_sample_barcodes(self):
        return self.sample_barcodes



# class DrugResponse():
#     def __init__(self, drugs_file_path="nationwidechildrens.org_clinical_drug.txt", patient_column="bcr_patient_barcode",
#                  columns=['bcr_patient_barcode', 'pharmaceutical_therapy_drug_name', 'pharmaceutical_therapy_type', 'treatment_best_response'],
#                  drug_name_col=None, response_column=None):
#         self.drug_name_col = drug_name_col
#         self.response_column = response_column
#
#         self.drugs = pd.read_table(drugs_file_path,
#                                    sep="\t",
#                                    skiprows=[1, 2],
#                                    na_values=["[Not Available]", "[Unknown]", "[Not Applicable]"],
#                                    usecols=columns
#                                    )
#         self.drugs.set_index(patient_column, inplace=True)


# class Biospecimen():
#     def __init__(self, biospecimens_file="genome.wustl.edu_biospecimen_sample.txt", patient_col_name="bcr_patient_barcode",
#                  columns=['bcr_sample_barcode', 'sample_type']):
#         self.biospecimen = pd.read_table(biospecimens_file, sep="\t", skiprows=[1, ],
#                                          na_values=["[Not Available]", "[Unknown]", "[Not Applicable]"],
#                                          usecols=columns
#                                          )
#         self.sample_barcodes = self.biospecimen[patient_col_name].tolist()
#         self.biospecimen.set_index(patient_col_name, inplace=True)
