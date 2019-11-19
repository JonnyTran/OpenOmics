import dash
import dash_bootstrap_components as dbc

from openomics_web.layouts import app_layout
from openomics_web.server import server

# running directly with Python
app = dash.Dash(__name__,
                server=server,
                external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = app_layout.app_main()

if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
