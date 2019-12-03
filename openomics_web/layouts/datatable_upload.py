import base64

import dash_html_components as html
import dash_table as dt

from openomics import MicroRNA, MessengerRNA


def parse_datatable_file(cohort_name, contents, filename, data_type):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        # with open('mytsvfile.tsv', 'r') as tsv:
        #     columns = tsv.readline().split('\t')
        if data_type == MicroRNA.name():
            df = MicroRNA(filename, decoded.decode('utf-8'))
        elif data_type == MessengerRNA.name():
            df = MessengerRNA(filename, decoded.decode('utf-8'))

        # # df = ExpressionData(filename, decoded.decode('utf-8'))
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
        ExpressionDataTable(df.head()),
    ])


def ExpressionDataTable(df):
    return dt.DataTable(
        id='table',
        columns=[{"name": i, "id": i} for i in df.columns],
        data=df.to_dict('records'),
        style_as_list_view=True,
        style_cell={'textAlign': 'left',
                    "maxWidth": 100, },
        style_data_conditional=[
            {'if': {'row_index': 'odd'},
             'backgroundColor': 'rgb(248, 248, 248)'
             },
        ],
        style_table={"maxHeight": "310px",
                     'width': '320px',
                     'marginTop': '5px',
                     'marginBottom': '10px',
                     },
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        row_selectable="multi",
        selected_rows=[],
        page_action="native",
        page_current=0,
        page_size=10,
    )


def expression_data_view():
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
