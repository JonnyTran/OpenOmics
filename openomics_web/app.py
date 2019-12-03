import dash
from dash.dependencies import Input, Output, State

from openomics_web.layouts import app_layout
from openomics_web.layouts.datatable_upload import parse_datatable_file
from openomics_web.server import server

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

# running directly with Python
app = dash.Dash(__name__,
                server=server,
                external_stylesheets=external_stylesheets)

app.layout = app_layout.app_main()


@app.callback([Output('data_table_columns_select', 'options')],
              [Input('upload-data', 'contents'),
               Input('upload-data', 'filename')],
              [State('data-table-type', 'value'), ])
def update_datatable_upload(list_of_contents, list_of_names, cohort_name, data_type, ):
    print("cohort_name", cohort_name)
    print("list_of_names", list_of_names)
    print("data_type", data_type)
    print("list_of_contents", list_of_contents)

    children = [parse_datatable_file(cohort_name, c, n, data_type) for c, n in zip(list_of_contents, list_of_names)]

    return


@app.callback(Output('output-data-upload', 'children'),
              [Input('submit-button', 'n_clicks')],
              [State('data-table-cohort', 'value'),
               State('data-table-type', 'value'),
               State('upload-data', 'contents'),
               State('upload-data', 'filename'),
               State('data_table_columns_select', 'value')
               ])
def import_datatable_upload(n_clicks, cohort_name, data_type, list_of_contents, list_of_names, columns):
    if list_of_contents is None:
        raise Exception()

    print("columns", columns)

    children = [parse_datatable_file(cohort_name, c, n, data_type) for c, n in zip(list_of_contents, list_of_names)]
    return children

if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
