import dash_core_components as dcc
import dash_html_components as html

from openomics_web.layouts.datatable_upload import datatable_upload


def control_tabs():
    return html.Div(id='circos-control-tabs', className='control-tabs', children=[
        dcc.Tabs(id='circos-tabs', value='what-is', children=[
            dcc.Tab(
                label='About',
                value='what-is',
                children=html.Div(className='control-tab', children=[
                    html.H4(className='what-is', children="What is Circos?"),

                    html.P('In the "Table" tab, you can view the datasets that define '
                           'the parameters of the graph, such as the layout, the '
                           'highlights, and the chords. You can interact with Circos '
                           'through this table by selecting the "Chords" graph in the '
                           '"Graph" tab, then viewing the "Chords" dataset in the '
                           '"Table" tab.'),

                    html.Div([
                        'Reference: ',
                        html.A('Seminal paper',
                               href='http://www.doi.org/10.1101/gr.092759.109)')
                    ]),
                    html.Div([
                        'For a look into Circos and the Circos API, please visit the '
                        'original repository ',
                        html.A('here', href='https://github.com/nicgirault/circosJS)'),
                        '.'
                    ]),

                    html.Br()
                ])
            ),

            dcc.Tab(
                label='Data',
                value='data',
                children=[datatable_upload()]
            )

        ])
    ])
