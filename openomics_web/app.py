import base64
import io
import os

import dash
import dash_html_components as html
import pandas as pd
from dash.dependencies import Input, Output, State

from openomics import MicroRNA, MessengerRNA, LncRNA, ProteinExpression
from openomics_web.layouts import app_layout
from openomics_web.layouts.datatable_upload import ExpressionDataTable, DataTableColumnSelect
from openomics_web.server import server

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

# running directly with Python
app = dash.Dash(__name__,
                server=server,
                external_stylesheets=external_stylesheets)

app.layout = app_layout.app_main()


@app.callback([Output('upload_table_preview', 'children'),
               Output('upload-data', 'children')],
              [Input('upload-data', 'contents'),
               Input('upload-data', 'filename')],
              [State('data-table-type', 'value'), ])
def update_datatable_metadata(list_of_contents, list_of_names, data_type, ):
    if list_of_contents is None:
        return None

    try:
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

    except Exception as e:
        print(e)
        return 'There was an error processing this file.'

    return DataTableColumnSelect(columns), "Uploaded {}".format(list_of_names)


@app.callback(Output('output-data-upload', 'children'),
              [Input('submit-button', 'n_clicks')],
              [State('data-table-cohort', 'value'),
               State('data-table-type', 'value'),
               State('upload-data', 'contents'),
               State('upload-data', 'filename'),
               State('data_table_genes_col_name', 'value'),
               State('data_table_columns_select', 'value')
               ])
def import_datatable_upload(n_clicks, cohort_name, data_type, list_of_contents, list_of_names, genes_col_name,
                            columns_select):
    if list_of_contents is None: return []

    columns = "|".join(columns_select)

    if len(list_of_contents) == 1:
        content = list_of_contents[0]
        filename = list_of_names[0]
    else:
        print("Multiple Files not supported")

    content_type, content_string = content.split(',')
    decoded = base64.b64decode(content_string)

    try:
        if 'csv' in filename:  # Assume that the user uploaded a CSV file
            file = io.StringIO(decoded.decode('utf-8'))
        elif 'xls' in filename:  # Assume that the user uploaded an excel file
            file = io.BytesIO(decoded)
        elif 'tsv' in filename:  # Assume that the user uploaded an tsv file
            file = io.StringIO(decoded.decode('utf-8'))
        elif 'txt' in filename:  # Assume that the user uploaded either a tsv or csv file
            file = decoded.decode('utf-8')

        if data_type == MicroRNA.name():
            df = MicroRNA(cohort_name, file, columns=columns, genes_col_name=genes_col_name)
        elif data_type == MessengerRNA.name():
            df = MessengerRNA(cohort_name, file, columns=columns, genes_col_name=genes_col_name)
        elif data_type == LncRNA.name():
            df = LncRNA(cohort_name, file, columns=columns, genes_col_name=genes_col_name)
        elif data_type == ProteinExpression.name():
            df = ProteinExpression(cohort_name, file, columns=columns, genes_col_name=genes_col_name)

    except Exception as e:
        # print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    return ExpressionDataTable(df)

if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
