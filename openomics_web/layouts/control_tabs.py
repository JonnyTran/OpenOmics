import dash_core_components as dcc
import dash_html_components as html

from openomics import MicroRNA, LncRNA, MessengerRNA, ProteinExpression


def control_tabs():
    return html.Div(id='circos-control-tabs', className='control-tabs', children=[
        dcc.Tabs(id='circos-tabs', value='what-is', children=[
            dcc.Tab(
                label='Welcome',
                value='what-is',
                children=html.Div(className='control-tab', children=[
                    html.H5(className='what-is', children="OpenOmics: Multi-Omics Data Explorer"),

                    html.P('OpenOmics provides a bioinformatics API and web-app'
                           'platform integrate, analyze, and visualize the'
                           'multi-omics and clinical data. You can interact with the'
                           'data input tabs to integrate the various different data types'
                           'including bio-databases, multi-omics expression, genomics,'
                           'and clinical data.'),

                    html.Div([
                        'Cite us at: ',
                        html.A('OpenOmics paper',
                               href='http://www.doi.org/10.1101/gr.092759.109'),
                        '.'
                    ]),
                    html.Div([
                        'For a look into the OpenOmics API and usage documentation, please visit the '
                        'original repository ',
                        html.A('OpenOmics GitHub', href='https://github.com/JonnyTran/OpenOmics.git'),
                        '.'
                    ]),

                    html.Br(),
                    html.Button(id='start-button', style={"align": "center"}, n_clicks=0,
                                children='Get Started'),
                ])
            ),

            dcc.Tab(
                label='Import Data Table',
                value='data',
                className='control-tab',
                children=[
                    html.Div(['Cohort name: ', ]),
                    dcc.Input(id='data-table-cohort', style={'width': '100%'},
                              placeholder="Leave empty to use filename",
                              type='text'),
                    html.Br(),
                    html.Div(['Select data type:', ]),
                    dcc.RadioItems(id='data-table-type',
                                   options=[
                                       {'label': 'Protein Expression', 'value': ProteinExpression.name()},
                                       {'label': 'miRNA Expression', 'value': MicroRNA.name()},
                                       {'label': 'lncRNA Expression', 'value': LncRNA.name()},
                                       {'label': 'mRNA Expression', 'value': MessengerRNA.name()},
                                   ],
                                   value=MicroRNA.name(), ),
                    html.Br(),
                    html.Div(['Import Data Table', ]),
                    html.Div(children=[
                        dcc.Upload(
                            id='upload-data',
                            children=html.Div(['Drag and Drop or ', html.A('Select Files')]),
                            style={
                                'width': '100%',
                                'height': '60px',
                                'lineHeight': '60px',
                                'borderWidth': '1px',
                                'borderStyle': 'dashed',
                                'borderRadius': '5px',
                                'textAlign': 'center',
                                'margin': '10px'
                            },
                            multiple=True
                        )
                    ]),
                    html.Br(),
                    html.Div(id='upload_table_preview'),
                    html.Br(),
                    html.Button(id='submit-button', n_clicks=0, children='Submit'),
                ]
            )

        ])
    ])
