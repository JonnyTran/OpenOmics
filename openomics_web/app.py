import dash
import dash_html_components as html
from dash.dependencies import Input, Output, State

from openomics_web.layouts import app_layout
from openomics_web.layouts.datatable_view import ExpressionDataTable, DataTableColumnSelect
from openomics_web.server import server
from openomics_web.utils.io import get_table_columns, get_expression_table

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
    if list_of_contents is None: return None
    try:
        columns = get_table_columns(list_of_contents, list_of_names)

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
               State('data-table-genes-col-name', 'value'),
               State('data-table-columns-select', 'value'),
               State('data-table-transpose', 'value')
               ])
def import_datatable_upload(n_clicks, cohort_name, data_type, list_of_contents, list_of_names, genes_col_name,
                            columns_select, transposed):
    if list_of_contents is None: return []
    try:
        df = get_expression_table(list_of_contents, list_of_names, data_type, cohort_name, genes_col_name,
                                  columns_select, transposed)
        print(df.describe())
    except Exception as e:
        print(e)
        return html.Div(['There was an error processing this file.'])

    return ExpressionDataTable(df)


if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
