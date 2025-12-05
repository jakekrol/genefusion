#!/usr/bin/env python3
import dash
from dash import dash_table, html
import pandas as pd
import argparse
parser = argparse.ArgumentParser(description="Display big data in a Dash DataTable")
parser.add_argument("-i", "--input", default="your_big_data.tsv", help="Input tsv")
parser.add_argument("-p", "--port", type=int, default=9999, help="Port to run the Dash app on")
args = parser.parse_args()

df = pd.read_csv(args.input, sep="\t")

app = dash.Dash(__name__)

app.layout = html.Div([
    dash_table.DataTable(
        data=df.head(500).to_dict("records"),
        columns=[{"name": c, "id": c} for c in df.columns],
        filter_action="native",
        sort_action="native",
        page_size=20
    )
])

if __name__ == "__main__":
    app.run(debug=True, port=args.port)
