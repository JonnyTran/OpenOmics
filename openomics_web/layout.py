import os

import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt


def app_layout():
    return html.Div(id='circos-body', className='app-body', children=[
        loading(),
        control_tabs(),
        table_view()
    ])


def loading():
    return html.Div(id='loading', children=[
        dcc.Loading(className='loading', children=html.Div(id="circos-hold", ))]
                    )


def table_view():
    return html.Div(id='table-container', children=[dt.DataTable(
        id="data-table",
        row_selectable='multi',
        # sorting=True,
        # filtering=True,
        css=[{
            "selector": ".dash-cell div.dash-cell-value",
            "rule": "display: inline; "
                    "white-space: inherit; "
                    "overflow: auto; "
                    "text-overflow: inherit;"
        }],
        style_cell={
            "whiteSpace": "no-wrap",
            "overflow": "hidden",
            "textOverflow": "ellipsis",
            "maxWidth": 100,
            'fontWeight': 100,
            'fontSize': '11pt',
            'fontFamily': 'Courier New',
            'backgroundColor': '#1F2132'
        },
        style_header={
            'backgroundColor': '#1F2132',
            'textAlign': 'center'
        },
        style_table={
            "maxHeight": "310px",
            'width': '320px',
            'marginTop': '5px',
            'marginBottom': '10px',
        },
        # n_fixed_rows=1,
        # n_fixed_columns=1
    )])


def control_tabs():
    return html.Div(id='circos-control-tabs', className='control-tabs', children=[
        dcc.Tabs(id='circos-tabs', value='what-is', children=[
            dcc.Tab(
                label='About',
                value='what-is',
                children=html.Div(className='control-tab', children=[
                    html.H4(className='what-is', children="What is Circos?"),

                    html.P('A Dash Circos graph consists of two main parts: the layout '
                           'and the tracks. '
                           'The layout sets the basic parameters of the graph, such as '
                           'radius, ticks, labels, etc; the tracks are graph layouts '
                           'that take in a series of data points to display.'),
                    html.P('The visualizations supported by Dash Circos are: heatmaps, '
                           'chords, highlights, histograms, line, scatter, stack, '
                           'and text graphs.'),
                    html.P('In the "Data" tab, you can opt to use preloaded datasets; '
                           'additionally, you can download sample data that you would '
                           'use with a Dash Circos component, upload that sample data, '
                           'and render it with the "Render" button.'),
                    html.P('In the "Graph" tab, you can choose the type of Circos graph '
                           'to display, control the size of the graph, and access data '
                           'that are generated upon hovering over parts of the graph. '),
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
                children=html.Div(className='control-tab', children=[
                    html.Div(className='app-controls-block', children=[
                        html.Div(className='app-controls-name', children='Data source'),
                        dcc.Dropdown(
                            id='circos-preloaded-uploaded',
                            options=[
                                {'label': 'Preloaded', 'value': 'preloaded'},
                                {'label': 'Upload', 'value': 'upload'}
                            ],
                            value='preloaded'
                        )
                    ]),
                    html.Hr(),
                    html.A(
                        html.Button(
                            id='circos-download-button',
                            className='control-download',
                            children="Download sample data"
                        ),
                        href=os.path.join('assets', 'sample_data', 'circos_sample_data.rar'),
                        download="circos_sample_data.rar",
                    ),

                    html.Div(id='circos-uploaded-data', children=[
                        dcc.Upload(
                            id="upload-data",
                            className='control-upload',
                            children=html.Div(
                                [
                                    "Drag and Drop or "
                                    "click to import "
                                    ".CSV file here!"
                                ]
                            ),
                            multiple=True,
                        ),
                        html.Div(className='app-controls-block', children=[
                            html.Div(className='app-controls-name',
                                     children='Select upload data'),
                            dcc.Dropdown(
                                id="circos-view-dataset-custom",
                                options=[
                                    {
                                        "label": "Layout",
                                        "value": 0,
                                    },
                                    {
                                        "label": "Track 1",
                                        "value": 1,
                                    },
                                    {
                                        "label": "Track 2",
                                        "value": 2,
                                    },
                                ],
                                value=0,
                            ),
                        ]),
                        html.Button(
                            "Render uploaded dataset",
                            id="render-button",
                            className='control-download',
                        )

                    ]),
                ])
            ),

            dcc.Tab(
                label='Graph',
                value='graph',
                children=html.Div(className='control-tab', children=[
                    html.Div(className='app-controls-block', children=[
                        html.Div(className='app-controls-name', children='Graph type'),
                        dcc.Dropdown(
                            id='circos-graph-type',
                            options=[
                                {'label': graph_type.title(),
                                 'value': graph_type} for graph_type in [
                                    'heatmap',
                                    'chords',
                                    'highlight',
                                    'histogram',
                                    'line',
                                    'scatter',
                                    'stack',
                                    'text',
                                    'parser_data'
                                ]
                            ],
                            value='chords'
                        ),
                        html.Div(className='app-controls-desc', id='chords-text'),
                    ]),
                    html.Div(className='app-controls-block', children=[
                        html.Div(className='app-controls-name', children='Graph size'),
                        dcc.Slider(
                            id='circos-size',
                            min=500,
                            max=800,
                            step=10,
                            value=650
                        ),
                    ]),
                    html.Hr(),
                    html.H5('Hover data'),
                    html.Div(
                        id='event-data-select'
                    ),

                ]),
            ),

        ])
    ])
