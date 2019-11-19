import dash_bio
import dash_core_components as dcc
import dash_html_components as html

from openomics_web.layouts.sidebar import control_tabs


def app_main():
    return html.Div(id='circos-body', className='app-body', children=[
        loading(),
        control_tabs(),
        html.Div(
            [
                html.Div(id="output-data-upload"),
                html.Div(id="event-data-store"),
            ],
            className="circos-display-none",
        ),
    ])



def loading():
    return html.Div(id='loading', children=[
        dcc.Loading(className='loading', children=html.Div(
            id="circos-hold",
            # children=[empty]
        ))
    ])


empty = dash_bio.Circos(
    id='main-circos',
    selectEvent={},
    layout=[],
    # size=700,
    config={},
    tracks=[],
    enableZoomPan=True,
    enableDownloadSVG=False
)
