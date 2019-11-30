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


@app.callback(Output('output-data-upload', 'children'),
              [Input('upload-data', 'contents'),
               Input('data-table-type', 'value')],
              [State('upload-data', 'filename'),
               State('upload-data', 'last_modified')])
def update_datatable_upload(list_of_contents, data_type, list_of_names, list_of_dates):
    if list_of_contents is not None:
        children = [
            parse_datatable_file(c, n, data_type) for c, n in
            zip(list_of_contents, list_of_names)]
        return children

if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
