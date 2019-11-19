import dash_html_components as html
import dash_table as dt


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
