import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt
import pandas as pd


def ClinicalDataColumnSelect(columns):
    return html.Div([
        html.Div(['Select the gene id/name column to index by:']),
        dcc.Dropdown(
            id='clinical-patient-col-name',
            options=[{'label': col, 'value': col} for col in columns],
            style={
                'width': '100%',
            },
            value=columns[0],
        ),
        html.Div(['Select the column prefixes to import:']),
        dcc.Dropdown(
            id='clinical-data-columns-select',
            options=[{'label': col, 'value': col} for col in columns],
            style={
                'width': '100%',
            },
            value=columns,
            multi=True,
        )
    ])


def ClinicalDataTable(df: pd.DataFrame):
    df.index.rename("id", inplace=True)

    return dt.DataTable(
        id='clinical-datatable',
        columns=[{"name": i, "id": i} for i in df.columns],
        data=df.reset_index().to_dict('records'),
        style_as_list_view=True,
        style_cell={'textAlign': 'left',
                    "maxWidth": 100, },
        style_data_conditional=[
            {'if': {'row_index': 'odd'},
             'backgroundColor': 'rgb(248, 248, 248)'
             },
        ],
        style_table={"maxHeight": '800px',
                     'width': '800px',
                     'marginTop': '5px',
                     'marginBottom': '10px',
                     'overflowX': 'scroll'
                     },
        virtualization=True,
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        row_selectable="multi",
        selected_rows=[],
        page_action="native",
        page_current=0,
        page_size=10,
    )
