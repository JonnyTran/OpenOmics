import base64
import io
import os

import dash
import pandas as pd
from dash.dependencies import Input, Output, State

from openomics_web.layouts import app_layout
from openomics_web.layouts.datatable_upload import parse_datatable_file, DataTableColumnSelect
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
    if list_of_contents is None:
        return []
    columns = "|".join(columns_select)
    print("columns", columns)
    print("genes_col_name", genes_col_name)
    return parse_datatable_file(cohort_name, list_of_contents, list_of_names, data_type, columns, genes_col_name)

if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
