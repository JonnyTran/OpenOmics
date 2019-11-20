import base64
import datetime

import dash_core_components as dcc
import dash_html_components as html

from openomics import ExpressionData
from openomics_web.views.table import ExpressionDataView


def datatable_upload():
    return html.Div(className='control-tab', children=[
        dcc.Upload(
            id='upload-data',
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select Files')
            ]),
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
            # Allow multiple files to be uploaded
            multiple=True
        ),

    ])


def parse_datatable_file(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        df = ExpressionData(filename, decoded.decode('utf-8'))
        # if 'csv' in filename:
        #     # Assume that the user uploaded a CSV file
        #     df = pd.read_csv(
        #         io.StringIO(decoded.decode('utf-8')))
        # elif 'xls' in filename:
        #     # Assume that the user uploaded an excel file
        #     df = pd.read_excel(io.BytesIO(decoded))
        # elif 'tsv' in filename:
        #     # Assume that the user uploaded an tsv file
        #     df = pd.read_table(io.StringIO(decoded.decode('utf-8')))
        # elif 'txt' in filename:
        #     # Assume that the user uploaded either a tsv or csv file
        #     df = pd.read_table(
        #         io.StringIO(decoded.decode('utf-8'))
        #     )

    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    return html.Div([
        html.H5(filename),
        html.H6(datetime.datetime.fromtimestamp(date)),

        ExpressionDataView(df.head()),

        html.Hr(),  # horizontal line

        # For debugging, display the raw contents provided by the web browser
        html.Div('Raw Content'),
        html.Pre(contents[0:200] + '...', style={
            'whiteSpace': 'pre-wrap',
            'wordBreak': 'break-all'
        })
    ])
