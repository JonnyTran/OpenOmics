import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt

from openomics_web.utils.str_utils import longest_common_prefix


def DataTableColumnSelect(columns):
    longest_common_prefixes = longest_common_prefix(columns)

    return html.Div([
        html.Div(['Select the gene id/name column to index by:']),
        dcc.Dropdown(
            id='data-table-genes-col-name',
            options=[{'label': col, 'value': col} for col in columns],
            style={
                'width': '100%',
            },
            value=columns[0],
        ),
        html.Div(['Select the column prefixes to import:']),
        dcc.Dropdown(
            id='data-table-columns-select',
            options=[{'label': col, 'value': col} for col in longest_common_prefixes],
            style={
                'width': '100%',
            },
            multi=True,
        )
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
