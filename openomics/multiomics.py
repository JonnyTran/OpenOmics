import logging
import os
import warnings
from os.path import exists, join
from typing import List, Dict, Union

import pandas as pd
from logzero import logger
from ruamel import yaml

import openomics
from .clinical import (
    ClinicalData,
    HISTOLOGIC_SUBTYPE_COL,
    PATHOLOGIC_STAGE_COL,
    TUMOR_NORMAL_COL,
    PREDICTED_SUBTYPE_COL,
)
from .database.base import Annotatable
from .genomics import SomaticMutation, CopyNumberVariation, DNAMethylation
from .imageomics import WholeSlideImage
from .proteomics import Protein
from .transcriptomics import MessengerRNA, MicroRNA, LncRNA, Expression

__all__ = ['MultiOmics']

class MultiOmics:
    """A data object which holds multiple -omics data for a single clinical
    cohort.
    """
    def __init__(self, cohort_name, omics_data=None):
        """
        Args:
            cohort_name (str): the clinical cohort name
            omics_data:
        """
        self._cohort_name = cohort_name
        self._omics = []

        # This is a data dictionary accessor to retrieve individual -omic data
        self.data: Dict[str, pd.DataFrame] = {}

        if omics_data:
            for omics in omics_data:
                self.add_omic(omics)

    def __repr__(self):
        return f"Cohort: {self._cohort_name}" \
               "\nExpression: {}" \
               "\nAnnotations: {}".format(
            {ntype: df.shape for ntype, df in self.data.items()},
            {ntype: omic.annotations.shape for ntype, omic in self.__dict__.items() if hasattr(omic, 'annotations')})

    @classmethod
    def load(cls, path):
        """
        Load a MultiOmics save directory and create a MultiOmics object.
        Args:
            path (str):

        Returns:
            MultiOmics
        """
        if isinstance(path, str) and '~' in path:
            path = os.path.expanduser(path)

        with open(join(path, 'metadata.yml'), 'r') as f:
            warnings.simplefilter('ignore', yaml.error.UnsafeLoaderWarning)
            attrs: Dict = yaml.load(f)

        logger.info(f"Loading MultiOmics {attrs} from {path}")

        omics = []
        for name in attrs['_omics']:
            if exists(join(path, f'{name}.expressions.pickle')):
                df = pd.read_pickle(join(path, f'{name}.expressions.pickle'))
                if name == MessengerRNA.name():
                    OmicsClass = MessengerRNA
                elif name == MicroRNA.name():
                    OmicsClass = MicroRNA
                elif name == LncRNA.name():
                    OmicsClass = LncRNA
                elif name == Protein.name():
                    OmicsClass = Protein
                elif name == Expression.name():
                    OmicsClass = Expression

                omic = OmicsClass(df, transpose=False)
            else:
                continue

            # annotation data
            if exists(join(path, f'{name}.annotations.pickle')):
                omic.annotations = pd.read_pickle(join(path, f'{name}.annotations.pickle'))

            omics.append(omic)

        self = cls(cohort_name=attrs['_cohort_name'], omics_data=omics)
        return self

    def save(self, path):
        """
        Serialize all omics data to disk.
        Args:
            path (str): A path to a directory where the MultiOmics data will be serialized into files.
        """
        if isinstance(path, str) and '~' in path:
            path = os.path.expanduser(path)

        if not exists(path):
            os.makedirs(path)

        # Write metadata to YAML so can be readable
        attrs = {
            '_omics': getattr(self, '_omics', None),
            '_cohort_name': getattr(self, '_cohort_name', None),
        }
        with open(join(path, 'metadata.yml'), 'w') as outfile:
            yaml.dump(attrs, outfile, default_flow_style=False)

        # Write each omics type
        for name, expressions_df in self.data.items():
            # expressions data
            if not exists(join(path, f'{name}.expressions.pickle')):
                expressions_df.to_pickle(join(path, f'{name}.expressions.pickle'))

            # annotation data
            expression_obj = self.__getitem__(name)
            if not exists(join(path, f'{name}.annotations.pickle')):
                expression_obj.annotations.to_pickle(join(path, f'{name}.annotations.pickle'))

    def add_omic(self,
                 omic_data: Union[Expression, Annotatable],
                 init_annotations: bool = True):
        """Adds an omic object to the Multiomics such that the samples in omic
        matches the samples existing in the other omics.

        Args:
            omic_data (Expression): The omic to add, e.g., MessengerRNA,
                MicroRNA, LncRNA, etc.
            init_annotations (bool): default True. If true, initializes
                the annotation dataframe in the omic object
        """
        self.__dict__[omic_data.name()] = omic_data

        if omic_data.name() not in self._omics:
            self._omics.append(omic_data.name())

        # dictionary as data accessor to the expression data
        self.data[omic_data.name()] = omic_data.expressions

        # Initialize annotation
        if init_annotations:
            omic_data.init_annotations(index=omic_data.get_genes_list())

        logging.info(
            "{} {} , indexed by: {}".format(omic_data.name(),
                                            self.data[omic_data.name()].shape if hasattr(
                                                self.data[omic_data.name()], "shape") else ": None",
                                            omic_data.annotations.index.name))

    def add_clinical_data(self, clinical: openomics.clinical.ClinicalData, **kwargs):
        """Add a ClinicalData instance to the MultiOmics instance.

        Args:
            clinical (openomics.clinical.ClinicalData):
            **kwargs:
        """
        if not isinstance(clinical, ClinicalData):
            raise Exception("Must pass a ClinicalData in, not a file path.")

        self.clinical = clinical

        self.data["PATIENTS"] = self.clinical.patient
        if hasattr(self.clinical, "biospecimen"):
            self.data["BIOSPECIMENS"] = self.clinical.biospecimen
        if hasattr(self.clinical, "drugs"):
            self.data["DRUGS"] = self.clinical.drugs

        self.build_samples(**kwargs)

    def get_omics_list(self):
        return self._omics

    def __getitem__(self, item: str) -> Union[Expression, Annotatable]:
        """This function allows the MultiOmicData class objects to access
        individual omics by a dictionary lookup, e.g. openomics["MicroRNA"]

        Args:
            item (str): a string of the class name
        """
        if item.lower() == MessengerRNA.name().lower():
            return self.__getattribute__(MessengerRNA.name())

        elif item.lower() == MicroRNA.name().lower():
            return self.__getattribute__(MicroRNA.name())

        elif item.lower() == LncRNA.name().lower():
            return self.__getattribute__(LncRNA.name())

        elif item.lower() == WholeSlideImage.name().lower():
            return self.__getattribute__(WholeSlideImage.name())

        elif item.lower() == SomaticMutation.name().lower():
            return self.__getattribute__(SomaticMutation.name())

        elif item.lower() == CopyNumberVariation.name().lower():
            return self.__getattribute__(CopyNumberVariation.name())

        elif item.lower() == DNAMethylation.name().lower():
            return self.__getattribute__(DNAMethylation.name())

        elif item.lower() == Protein.name().lower():
            return self.__getattribute__(Protein.name())

        elif item.lower() == "patients":
            return self.clinical.patient
        elif item.lower() == "samples":
            if hasattr(self, "clinical"):
                return self.clinical.samples
            else:
                return self.samples
        elif item.lower() == "drugs":
            return self.clinical.drugs
        else:
            raise Exception(
                'String accessor must be one of {"MessengerRNA", "MicroRNA", "LncRNA", "Protein", etc.}'
            )

    def remove_duplicate_genes(self):
        """Removes duplicate genes between any omics such that the gene index
        across all omics has no duplicates.
        """
        for omic_A in self._omics:
            for omic_B in self._omics:
                if omic_A != omic_B:
                    self.__getattribute__(omic_A).drop_genes(
                        set(self.__getattribute__(omic_A).get_genes_list())
                        & set(self.__getattribute__(omic_B).get_genes_list()))

    def build_samples(self, agg_by="union"):
        """Running this function will build a dataframe for all samples across
        the different omics (either by a union or intersection). Then,

        Args:
            agg_by (str): ["union", "intersection"]
        """
        # make sure at least one ExpressionData present
        if len(self._omics) < 1:
            logging.debug(
                "build_samples() does nothing. Must add at least one omic to this MultiOmics object."
            )
            return

        all_samples = pd.Index([])
        for omic in self._omics:
            if agg_by == "union":
                all_samples = all_samples.union(self.data[omic].index)
            elif agg_by == "intersection":
                all_samples = all_samples.intersection(self.data[omic].index)

        if hasattr(self, "clinical"):
            self.clinical.build_clinical_samples(all_samples)
            self.samples = self.clinical.samples.index
        else:
            self.samples = all_samples

    def __dir__(self):
        return list(self.data.keys())

    def match_samples(self, omics) -> pd.Index:
        """Return the index of sample IDs of the intersection of
        samples from all modalities

        Args:
            omics: An array of modalities

        Returns:
            matched_sapmles: An pandas Index list
        """
        # TODO check that for single modalities, this fetch all patients
        matched_samples = self.data[omics[0]].index.copy()

        for omic in omics:
            matched_samples = matched_samples.join(self.data[omic].index,
                                                   how="inner")

        return matched_samples

    def load_data(self,
                  omics: Union[List[str], str],
                  target: List[str] = ["pathologic_stage"],
                  pathologic_stages=None,
                  histological_subtypes=None,
                  predicted_subtypes=None,
                  tumor_normal=None,
                  samples_barcode=None,
                  remove_duplicates=True):
        """ Prepare the multiomics data in format

        Args:
            omics (list): A list of the data modalities to load. Default "all"
                to select all modalities
            target (list): The clinical data fields to include in the
            pathologic_stages (list): Only fetch samples having certain stages
                in their corresponding patient's clinical data. For instance,
                ["Stage I", "Stage II"] will only fetch samples from Stage I and
                Stage II patients. Default is [] which fetches all pathologic
                stages.
            histological_subtypes: A list specifying the histological subtypes
                to fetch. Default is [] which fetches all histological sybtypes.
            predicted_subtypes: A list specifying the histological subtypes to
                fetch. Default is [] which fetches all histological sybtypes.
            tumor_normal: ["Tumor", "Normal"]. Default is [], which fetches all
                tumor or normal sample types.
            samples_barcode: A list of sample's barcode. If not None, only fetch
                data with matching samples provided in this list.
            remove_duplicates (bool): If True, only selects samples with non-duplicated index.

        Returns:
            Tuple[Dict[str, pd.DataFrame], pd.DataFrame]: Returns (X, y), where
            X is a dictionary containing the multiomics data with matched
            samples, and y contain the :param target: labels for those samples.
        """
        if omics == "all" or omics is None:
            omics = self._omics

        matched_samples = self.match_samples(omics)

        if samples_barcode is not None:
            matched_samples = samples_barcode

        if hasattr(self, "clinical") and isinstance(self.clinical,
                                                    ClinicalData):
            # Build targets clinical data
            y = self.get_sample_attributes(matched_samples)

            # Select only samples with certain cancer stage or subtype
            if pathologic_stages:
                y = y[y[PATHOLOGIC_STAGE_COL].isin(pathologic_stages)]
            if histological_subtypes:
                y = y[y[HISTOLOGIC_SUBTYPE_COL].isin(histological_subtypes)]
            if predicted_subtypes:
                y = y[y[PREDICTED_SUBTYPE_COL].isin(predicted_subtypes)]
            if tumor_normal:
                y = y[y[TUMOR_NORMAL_COL].isin(tumor_normal)]

            # Filter y target column labels
            y = y.filter(target)
            y.dropna(axis=0, inplace=True)
            matched_samples = y.index
        else:
            y = None

        # Build expression matrix for each omic, indexed by matched_samples
        X_multiomics = {}
        for omic in omics:
            X_multiomics[omic] = self.data[omic].loc[matched_samples, self[omic].get_genes_list()]

            if remove_duplicates:
                X_multiomics[omic] = X_multiomics[omic].loc[~X_multiomics[omic].index.duplicated(keep="first")]

        return X_multiomics, y

    def get_sample_attributes(self, matched_samples):
        """Fetch patient's clinical data for each given samples barcodes in the
        matched_samples

        Returns
            samples_index: Index of samples

        Args:
            matched_samples: A list of sample barcodes
        """
        return self.data["SAMPLES"].reindex(matched_samples)

    def print_sample_sizes(self):
        for omic in self.data:
            print(
                omic,
                self.data[omic].shape
                if hasattr(self.data[omic], "shape") else "None",
            )

    def annotate_samples(self, dictionary):
        """This function adds a "predicted_subtype" field to the patients
        clinical data. For instance, patients were classified into subtypes
        based on their expression profile using k-means, then, to use this
        function, do:

        annotate_patients(dict(zip(patient index>, <list of corresponding
        patient's subtypes>)))

        Adding a field to the patients clinical data allows openomics to
        query the patients data through the .load_data(subtypes=[]) parameter,

        Args:
            dictionary: A dictionary mapping patient's index to a subtype
        """
        self.data["PATIENTS"] = self.data["PATIENTS"].assign(
            subtypes=self.data["PATIENTS"][
                self.clinical.patient_column].map(dictionary))
