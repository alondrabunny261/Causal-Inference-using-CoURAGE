import marimo

__generated_with = "0.14.10-dev0"
app = marimo.App(width="medium")


@app.cell
def _():
    import pandas as pd
    import matplotlib.pyplot as plt
    return (pd,)


@app.cell
def _(pd):
    data = pd.read_csv("M1_30_min_meteorological_CO2_aerosol_data.csv")
    data
    return


if __name__ == "__main__":
    app.run()
