import base64
import io
import os

import pandas as pd

from openomics import MicroRNA, MessengerRNA, LncRNA, ProteinExpression


def get_table_columns(list_of_contents, list_of_names):
    content = list_of_contents[0]
    filename = list_of_names[0]
    _, file_extension = os.path.splitext(filename)
    content_type, content_string = content.split(',')
    decoded = base64.b64decode(content_string)
    first_line = decoded.decode('utf-8').partition('\n')[0]
    if file_extension == ".tsv":
        columns = first_line.split('\t')
    elif file_extension == ".csv":
        columns = first_line.split(',')
    elif file_extension == ".xls":
        columns = pd.read_excel(io.BytesIO(decoded), low_memory=True).columns.tolist()
    else:
        columns = pd.read_table(io.BytesIO(decoded), low_memory=True).columns.tolist()
    return columns


def get_expression_table(list_of_contents, list_of_names, data_type, cohort_name=None, genes_col_name=None,
                         columns_selected=None, transposed=None):
    """

    Args:
        list_of_contents:
        list_of_names:
        data_type:
        cohort_name:
        genes_col_name:
        columns_selected:

    Returns:

    """
    if columns_selected:
        columns = "|".join(columns_selected)
    else:
        columns = None

    if transposed == "True":
        transposed = True
    elif transposed == "False":
        transposed = False

    if len(list_of_contents) == 1:
        content = list_of_contents[0]
        filename = list_of_names[0]
    else:
        raise Exception("Multiple files not supported")  # TODO

    content_type, content_string = content.split(',')
    decoded = base64.b64decode(content_string)

    if 'csv' in filename:  # Assume that the user uploaded a CSV file
        file = decoded.decode('utf-8')
    elif 'xls' in filename:  # Assume that the user uploaded an excel file
        file = io.BytesIO(decoded)
    elif 'tsv' in filename:  # Assume that the user uploaded an tsv file
        file = decoded.decode('utf-8')
    elif 'txt' in filename:  # Assume that the user uploaded either a tsv or csv file
        file = decoded.decode('utf-8')

    if data_type == MicroRNA.name():
        df = MicroRNA(cohort_name, file, columns=columns, genes_col_name=genes_col_name, transposed=transposed)
    elif data_type == MessengerRNA.name():
        df = MessengerRNA(cohort_name, file, columns=columns, genes_col_name=genes_col_name, transposed=transposed)
    elif data_type == LncRNA.name():
        df = LncRNA(cohort_name, file, columns=columns, genes_col_name=genes_col_name, transposed=transposed)
    elif data_type == ProteinExpression.name():
        df = ProteinExpression(cohort_name, file, columns=columns, genes_col_name=genes_col_name, transposed=transposed)

    return df
