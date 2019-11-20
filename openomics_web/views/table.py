import dash_html_components as html
import dash_table as dt


def ExpressionDataView(df):
    return dt.DataTable(
        id='table',
        columns=[{"name": i, "id": i} for i in df.columns],
        data=df.to_dict('records'),
        style_as_list_view=True,
        style_cell={'textAlign': 'left'},
        style_data_conditional=[
            {
                'if': {'row_index': 'odd'},
                'backgroundColor': 'rgb(248, 248, 248)'
            }
        ],
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
