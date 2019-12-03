import base64
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
            columns = pd.read_excel(decoded.decode('utf-8')).columns.tolist()
        else:
            columns = first_line.split(',')

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
               # State('data_table_columns_select', 'value')
               ])
def import_datatable_upload(n_clicks, cohort_name, data_type, list_of_contents, list_of_names):
    if list_of_contents is None:
        return []

    # print("columns", columns)

    children = [parse_datatable_file(cohort_name, c, n, data_type) for c, n in zip(list_of_contents, list_of_names)]
    return children

if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
