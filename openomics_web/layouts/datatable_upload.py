import base64
import io

import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt

from openomics import MicroRNA, MessengerRNA, LncRNA, ProteinExpression
from openomics_web.utils.str_utils import longest_common_prefix


def DataTableColumnSelect(columns):
    longest_common_prefixes = longest_common_prefix(columns)
    print("longest_common_prefixes", longest_common_prefixes)

    return html.Div([
        html.Div(['Select the gene id/name column to index by:']),
        dcc.Dropdown(
            id='data_table_genes_col_name',
            options=[{'label': col, 'value': col} for col in columns],
            style={
                'width': '100%',
            },
        ),
        html.Div(['Select the column prefixes to import:']),
        dcc.Dropdown(
            id='data_table_columns_select',
            options=[{'label': col, 'value': col} for col in longest_common_prefixes],
            style={
                'width': '100%',
            },
            multi=True,
        )
    ])


def parse_datatable_file(cohort_name, contents, filenames, data_type, columns, genes_col):
    if len(contents) == 1:
        content = contents[0]
        filename = filenames[0]
    else:
        print("Multiple Files not supported")

    content_type, content_string = content.split(',')
    decoded = base64.b64decode(content_string)

    try:
        if 'csv' in filename:  # Assume that the user uploaded a CSV file
            file = io.StringIO(decoded.decode('utf-8'))
        elif 'xls' in filename:  # Assume that the user uploaded an excel file
            file = io.BytesIO(decoded)
        elif 'tsv' in filename:  # Assume that the user uploaded an tsv file
            file = io.StringIO(decoded.decode('utf-8'))
        elif 'txt' in filename:  # Assume that the user uploaded either a tsv or csv file
            file = decoded.decode('utf-8')

        if data_type == MicroRNA.name():
            df = MicroRNA(cohort_name, file, columns=columns, genes_col_name=genes_col)
        elif data_type == MessengerRNA.name():
            df = MessengerRNA(cohort_name, file, columns=columns, genes_col_name=genes_col)
        elif data_type == LncRNA.name():
            df = LncRNA(cohort_name, file, columns=columns, genes_col_name=genes_col)
        elif data_type == ProteinExpression.name():
            df = LncRNA(cohort_name, file, columns=columns, genes_col_name=genes_col)

        print("df.head()", df.head())

        # # df = ExpressionData(filename, decoded.decode('utf-8'))


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
