import marimo

__generated_with = "0.14.12-dev0"
app = marimo.App(width="full")


@app.cell
async def _():
    import pandas as pd
    import matplotlib.pyplot as plt
    from ipywidgets import interact, widgets
    import marimo as mo
    from IPython.display import display
    import plotly.express as px
    import micropip
    import json
    import numpy as np
    import datetime as dt
    import ast
    import plotly.graph_objects as go
    from scipy.stats import gaussian_kde

    await micropip.install("altair")
    import altair as alt
    return alt, ast, go, mo, pd, px


@app.cell
def _():
    ######################################################## Data sets ################################################################
    return


@app.cell
def _(pd):
    ### 30 Minute data ######
    M1_path = "M1_data"
    data_M1_30_min = pd.read_csv(
        M1_path
        + "/30_min_CO2_aerosols/M1_30_min_meteorological_CO2_aerosol_data.csv"
    )
    data_M1_30_min["specific humidity"] = data_M1_30_min[
        "co2 spec. humidity (meas)"
    ]
    data_M1_30_min["specific humidity"] = (
        data_M1_30_min["specific humidity"] * 1000
    )

    S2_path = "S2_data"
    data_S2_30_min = pd.read_csv(
        S2_path
        + "/30_min_CO2_aerosols/S2_30_min_meteorological_CO2_aerosol_data.csv"
    )
    data_S2_30_min["specific humidity"] = data_S2_30_min[
        "co2 spec. humidity (meas)"
    ]
    data_S2_30_min["specific humidity"] = (
        data_S2_30_min["specific humidity"] * 1000
    )

    S3_path = "S3_data"
    data_S3_30_min = pd.read_csv(
        S3_path
        + "/30_min_CO2_aerosols/S3_30_min_meteorological_CO2_aerosol_data.csv"
    )

    data_S3_30_min["specific humidity"] = data_S3_30_min[
        "co2 spec. humidity (meas)"
    ]
    data_S3_30_min["specific humidity"] = (
        data_S3_30_min["specific humidity"] * 1000
    )
    return (
        M1_path,
        S2_path,
        S3_path,
        data_M1_30_min,
        data_S2_30_min,
        data_S3_30_min,
    )


@app.cell
def _(M1_path, S2_path, S3_path, pd):
    ### 1 Minute data ######
    data_M1_1_min = pd.read_csv(
        M1_path
        + "/1_min_aerosols/M1_minute_meteorological_with_aerosols_data_wxt.csv"
    )

    data_S2_1_min = pd.read_csv(
        S2_path
        + "/1_min_aerosols/S2_minute_meteorological_with_aerosols_data_wxt.csv"
    )

    data_S3_1_min = pd.read_csv(
        S3_path
        + "/1_min_aerosols/S3_minute_meteorological_with_aerosols_data_wxt.csv"
    )
    return data_M1_1_min, data_S2_1_min, data_S3_1_min


@app.cell
def _(
    data_M1_1_min,
    data_M1_30_min,
    data_S2_1_min,
    data_S2_30_min,
    data_S3_1_min,
    data_S3_30_min,
):
    all_datasets = [
        data_M1_1_min,
        data_S2_1_min,
        data_S3_1_min,
        data_M1_30_min,
        data_S2_30_min,
        data_S3_30_min,
    ]
    return (all_datasets,)


@app.cell
def _(data_M1_30_min, data_S2_30_min, data_S3_30_min):
    min_30_datasets = [data_M1_30_min, data_S2_30_min, data_S3_30_min]
    return (min_30_datasets,)


@app.cell
def _(all_datasets, pd):
    for df in all_datasets:

        df["time"] = pd.to_datetime(df["time"])
        df["date"] = df["time"].dt.date
    return


@app.cell
def _(min_30_datasets):
    for element in min_30_datasets:

        element["total_N_conc"] = element["total number conc"]
    return


@app.cell
def _(pd):
    ################################### Morning #########################################
    M1_CI_coefs_30_min_m = pd.read_csv(
        "M1_data/30_min_CO2_aerosols/Morning/M1_significant_lags_30_min_C02_aerosols_morning_coefs.csv"
    )

    S2_CI_coefs_30_min_m = pd.read_csv(
        "S2_data/30_min_CO2_aerosols/Morning/S2_significant_lags_30_min_C02_aerosols_morning_coefs.csv"
    )

    S3_CI_coefs_30_min_m = pd.read_csv(
        "S3_data/30_min_CO2_aerosols/Morning/S3_significant_lags_30_min_C02_aerosols_morning_coefs.csv"
    )

    ################################### Night #########################################

    M1_CI_coefs_30_min_n = pd.read_csv(
        "M1_data/30_min_CO2_aerosols/Night/M1_significant_lags_30_min_C02_aerosols_night_coefs.csv"
    )

    S2_CI_coefs_30_min_n = pd.read_csv(
        "S2_data/30_min_CO2_aerosols/Night/S2_significant_lags_30_min_C02_aerosols_night_coefs.csv"
    )

    S3_CI_coefs_30_min_n = pd.read_csv(
        "S3_data/30_min_CO2_aerosols/Night/S3_significant_lags_30_min_C02_aerosols_night_coefs.csv"
    )

    ################################### Rain Events #########################################

    M1_CD_coefs_1_min = pd.read_csv(
        "M1_data/1_min_aerosols/M1_significant_lags_min_aerosols_coefs.csv"
    )

    S2_CD_coefs_1_min = pd.read_csv(
        "S2_data/1_min_aerosols/S2_significant_lags_min_aerosols_coefs.csv"
    )

    S3_CD_coefs_1_min = pd.read_csv(
        "S3_data/1_min_aerosols/S3_significant_lags_min_aerosols_coefs.csv"
    )
    return (
        M1_CD_coefs_1_min,
        M1_CI_coefs_30_min_m,
        M1_CI_coefs_30_min_n,
        S2_CD_coefs_1_min,
        S2_CI_coefs_30_min_m,
        S2_CI_coefs_30_min_n,
        S3_CD_coefs_1_min,
        S3_CI_coefs_30_min_m,
        S3_CI_coefs_30_min_n,
    )


@app.cell
def _():
    ############### Percentage of relation in events ###############
    return


@app.cell
def _(pd):
    ################################### Morning #########################################
    M1_CI_percentage_30_min_m = pd.read_csv(
        "M1_data/30_min_CO2_aerosols/Morning/M1_lag_summary_morning.csv"
    )

    S2_CI_percentage_30_min_m = pd.read_csv(
        "S2_data/30_min_CO2_aerosols/Morning/S2_lag_summary_morning.csv"
    )

    S3_CI_percentage_30_min_m = pd.read_csv(
        "S3_data/30_min_CO2_aerosols/Morning/S3_lag_summary_morning.csv"
    )

    ################################### Night #########################################
    M1_CI_percentage_30_min_m = pd.read_csv(
        "M1_data/30_min_CO2_aerosols/Night/M1_lag_summary_night.csv"
    )

    S2_CI_percentage_30_min_m = pd.read_csv(
        "S2_data/30_min_CO2_aerosols/Night/S2_lag_summary_night.csv"
    )

    S3_CI_percentage_30_min_m = pd.read_csv(
        "S3_data/30_min_CO2_aerosols/Night/S3_lag_summary_night.csv"
    )

    ################################### Rain #########################################
    M1_CI_percentage_1_min = pd.read_csv(
        "M1_data/1_min_aerosols/M1_lag_summary_rain.csv"
    )

    S2_CI_percentage_1_min = pd.read_csv(
        "S2_data/1_min_aerosols/S2_lag_summary_rain.csv"
    )

    S3_CI_percentage_1_min = pd.read_csv(
        "S3_data/1_min_aerosols/S3_lag_summary_rain.csv"
    )
    return (
        M1_CI_percentage_1_min,
        M1_CI_percentage_30_min_m,
        S2_CI_percentage_1_min,
        S2_CI_percentage_30_min_m,
        S3_CI_percentage_1_min,
        S3_CI_percentage_30_min_m,
    )


@app.cell
def _():
    ##################################################### Causal Estimation #########################################################
    return


@app.cell
def _(pd):
    ################################### Morning #########################################
    M1_WPC_30_min_m = pd.read_csv(
        "M1_data/30_min_CO2_aerosols/Morning/WPC_M1_morning_table.csv"
    )

    S2_WPC_30_min_m = pd.read_csv(
        "S2_data/30_min_CO2_aerosols/Morning/WPC_S2_morning_table.csv"
    )

    S3_WPC_30_min_m = pd.read_csv(
        "S3_data/30_min_CO2_aerosols/Morning/WPC_S3_morning_table.csv"
    )


    ################################### Night #########################################
    M1_WPC_30_min_n = pd.read_csv(
        "M1_data/30_min_CO2_aerosols/Night/WPC_M1_evening_table.csv"
    )

    S2_WPC_30_min_n = pd.read_csv(
        "S2_data/30_min_CO2_aerosols/Night/WPC_S2_night_table.csv"
    )

    S3_WPC_30_min_n = pd.read_csv(
        "S3_data/30_min_CO2_aerosols/Night/WPC_S3_evening_table.csv"
    )

    ################################### Rain #########################################
    M1_WPC_1_min = pd.read_csv(
        "M1_data/1_min_aerosols/WPC_M1_rain_event_table.csv"
    )

    S2_WPC_1_min = pd.read_csv(
        "S2_data/1_min_aerosols/WPC_S2_rain_event_table.csv"
    )

    S3_WPC_1_min = pd.read_csv(
        "S3_data/1_min_aerosols/WPC_S3_rain_event_table.csv"
    )
    return (
        M1_WPC_1_min,
        M1_WPC_30_min_m,
        M1_WPC_30_min_n,
        S2_WPC_1_min,
        S2_WPC_30_min_m,
        S2_WPC_30_min_n,
        S3_WPC_1_min,
        S3_WPC_30_min_m,
        S3_WPC_30_min_n,
    )


@app.cell
def _():
    Causal_Relations = [
        "Temperature-->Temperature",
        "Temperature-->Precipitation",
        "Temperature-->Specific Humidity",
        "Temperature-->Aerosols",
        "Temperature-->CO2",
        "Precipitation-->Temperature",
        "Precipitation-->Precipitation",
        "Precipitation-->Specific Humidity",
        "Precipitation-->Aerosols",
        "Precipitation-->CO2",
        "Specific Humidity-->Temperature",
        "Specific Humidity-->Precipitation",
        "Specific Humidity-->Specific Humidity",
        "Specific Humidity-->Aerosols",
        "Specific Humidity-->CO2",
        "Aerosols-->Temperature",
        "Aerosols-->Precipitation",
        "Aerosols-->Specific Humidity",
        "Aerosols-->Aerosols",
        "Aerosols-->CO2",
        "CO2-->Temperature",
        "CO2-->Precipitation",
        "CO2-->Specific Humidity",
        "CO2-->Aerosols",
        "CO2-->CO2",
    ]

    Causal_Relations_ = [
        "T_T",
        "T_P",
        "T_H",
        "T_aero",
        "T_CO2",
        "P_T",
        "P_P",
        "P_H",
        "P_aero",
        "P_CO2",
        "H_T",
        "H_P",
        "H_H",
        "H_aero",
        "H_CO2",
        "aero_T",
        "aero_P",
        "aero_H",
        "aero_aero",
        "aero_CO2",
        "CO2_T",
        "CO2_P",
        "CO2_H",
        "CO2_aero",
        "CO2_CO2",
    ]

    WPC_Relations_1_min = [
        "Precipitation on Aerosols",
        "Specific Humidity on Temperature",
        "Specific Humidity on Aerosols",
        "Specific Humidity on Precipitation",
        "Aerosols on Temperature",
    ]

    WPC_Relations_30_min = [
        "CO2 on Temperature",
        "Aerosols on Temperature",
        "Specific Humidity on Temperature",
        "Specific Humidity on Aerosols",
    ]
    return (
        Causal_Relations,
        Causal_Relations_,
        WPC_Relations_1_min,
        WPC_Relations_30_min,
    )


@app.cell
def _():
    return


@app.cell
def _():
    ############################### Functions #########################################
    return


@app.cell
def _(mo):
    def get_date():

        date_plot = mo.ui.date()

        return date_plot
    return (get_date,)


@app.cell
def _(get_date):
    date_plot = get_date()
    return (date_plot,)


@app.function
def select_date(data, value):
    selected_date = value
    filtered_data = data[data["time"].dt.date == selected_date]

    return filtered_data.head(1000)


@app.cell
def _(date_plot):
    selected_date = date_plot.value
    return (selected_date,)


@app.cell
def _(alt, mo):
    def plot_data(data, variable):
        if variable == "Precipitation":
            chart = mo.center(
                mo.ui.altair_chart(
                    alt.Chart(data, width=750)
                    .mark_circle()
                    .encode(
                        x=alt.X("time", axis=alt.Axis(labels=False)),
                        y="precipitation",
                    )
                )
            )
        elif variable == "Specific Humidity":
            chart = mo.center(
                mo.ui.altair_chart(
                    alt.Chart(data, width=750)
                    .mark_circle()
                    .encode(
                        x=alt.X("time", axis=alt.Axis(labels=False)),
                        y="specific humidity",
                    )
                )
            )

        elif variable == "Temperature":
            chart = mo.center(
                mo.ui.altair_chart(
                    alt.Chart(data, width=750)
                    .mark_circle()
                    .encode(
                        x=alt.X("time", axis=alt.Axis(labels=False)),
                        y="temperature",
                    )
                )
            )
        elif variable == "Aerosols":
            chart = mo.center(
                mo.ui.altair_chart(
                    alt.Chart(data, width=750)
                    .mark_circle()
                    .encode(
                        x=alt.X("time", axis=alt.Axis(labels=False)),
                        y="total_N_conc",
                    )
                )
            )

        return chart
    return (plot_data,)


@app.cell
def _(CI_histograms, Causal_Relations, Causal_Relations_):
    def select_the_CI_of_choice(sel_value, coef_table, facility, lag):

        if sel_value in Causal_Relations:
            position = Causal_Relations.index(sel_value)
            short_label = Causal_Relations_[position]
            coef_histogram = CI_histograms(
                coef_table, short_label, "magenta", short_label, facility, lag
            )

        else:
            coef_histogram = None

        return coef_histogram
    return (select_the_CI_of_choice,)


@app.cell
def _(Causal_Relations, Causal_Relations_, mo, pd, px):
    def CI_percentage_display(sel_value, data, lag, names):

        if sel_value in Causal_Relations:

            position = Causal_Relations.index(sel_value)

            short_label = Causal_Relations_[position].replace("_", " → ")

            if short_label in data.Relation.values:

                index = data.index[
                    (data["Relation"] == short_label) & (data["Lag"] == 1)
                ].tolist()

                if index:
                    i = index[0]  # ← primer valor del índice

                    sample_size = data["Total Samples"][i]

                    var_event_frequency = data["Count"][i]

                    percentage = data["Percentage"][i]

                    rest = sample_size - var_event_frequency

                    # Dataframe
                    dA = pd.DataFrame(
                        {
                            "Section": ["Sample", "Rest"],
                            "Value": [var_event_frequency, rest],
                        }
                    )

                    # Asegúrate de que no haya espacios u otras inconsistencias
                    dA["Section"] = dA["Section"].astype(str).str.strip()

                    # Define manualmente el orden de las categorías
                    dA["Section"] = pd.Categorical(
                        dA["Section"], categories=["Sample", "Rest"], ordered=True
                    )

                    # Donut graph
                    percentage_graph = px.pie(
                        dA,
                        values="Value",
                        names="Section",
                        hole=0.5,
                        color="Section",
                        color_discrete_map={
                            "Sample": "magenta",
                            "Rest": "#d3d3d3",
                        },
                        width=400,
                        height=400,
                    )

                    # Take away legend
                    percentage_graph.update_traces(textinfo="none")
                    percentage_graph.update_layout(showlegend=False)

                    # Update display
                    percentage_graph.update_layout(
                        annotations=[
                            dict(
                                text=f"{percentage}%",  # Texto personalizado
                                x=0.5,
                                y=0.5,
                                font_size=24,
                                showarrow=False,
                                font=dict(color="black"),
                            )
                        ]
                    )

                    simple_stat = mo.stat(
                        f"{var_event_frequency}/{sample_size}",
                        label="Causal Relation",
                        caption=f"of {sel_value}",
                        direction="increase",
                        bordered=True,
                    )

            else:

                percentage_graph = None

                simple_stat = mo.stat(
                    "0/100",
                    label="Causal Relation",
                    caption=f"of {sel_value}",
                    direction="increase",
                    bordered=True,
                )

        else:

            percentage_graph = None

            simple_stat = mo.stat(
                "0/100",
                label="Causal Relation",
                caption=f"of {sel_value}",
                direction="increase",
                bordered=True,
            )

        return simple_stat, percentage_graph
    return (CI_percentage_display,)


@app.cell
def _(ast, go, px):
    # Causal Estimation histograms
    def CE_histograms(data, selection, facility):

        # Aplica a toda la columna (reemplaza 'columna' con el nombre real)
        data[selection] = data[selection].apply(
            lambda x: ast.literal_eval(x)[0] if isinstance(x, str) else x
        )

        # Create interactive histograms
        fig = px.histogram(data[selection], width=700, height=500)
        mean_value = data[selection].mean()
        mean_str = "{:.3g}".format(mean_value)
        mean_value = float(mean_str)

        data[f"{selection} interval 1"] = data[f"{selection} interval 1"].apply(
            lambda x: ast.literal_eval(x)[0] if isinstance(x, str) else x
        )

        data[f"{selection} interval 2"] = data[f"{selection} interval 1"].apply(
            lambda x: ast.literal_eval(x)[0] if isinstance(x, str) else x
        )

        fig.add_trace(
            go.Histogram(
                x=data[f"{selection} interval 1"],
                name="Upper",
                marker_color="magenta",
            )
        )

        fig.add_trace(
            go.Histogram(
                x=data[f"{selection} interval 2"],
                name="Lower",
                marker_color="skyblue",
            )
        )

        # Mark median
        fig.add_vline(
            x=mean_value,
            line_dash="dash",
            line_color="red",
            annotation_text=f"Median = {mean_value}",
            annotation_position="top right",
        )

        # Add title
        fig.update_layout(
            title=facility,
            xaxis_title="Wright's path Estimation",
            yaxis_title="Frequency",
        )

        fig.show()

        return fig
    return (CE_histograms,)


@app.cell
def _(px):
    # Causal Discovery Histograms
    def CI_histograms(data_coef, combination, color, title, facility, lag):

        filtered = data_coef[
            (data_coef["var_pair"] == combination) & (data_coef["lag"] == lag)
        ]

        data = filtered["coef"]

        ##########################################################################
        # Create interactive histograms
        fig = px.histogram(
            data, nbins=20, color_discrete_sequence=[color], width=800, height=500
        )
        mean_value = data.mean()
        mean_str = "{:.3g}".format(mean_value)
        mean_value = float(mean_str)

        # Mark median
        fig.add_vline(
            x=mean_value,
            line_dash="dash",
            line_color="red",
            annotation_text=f"Median = {mean_value}",
            annotation_position="top right",
        )

        fig.add_vline(
            x=0,
            line_dash="dash",
            line_color="gray",
        )

        # Adjust margins
        fig.update_layout(
            title=f"{facility}: {title}",
            xaxis_range=[-1, 1],
            xaxis_title="Correlation Coefficients",
            yaxis_title="Frecuency",
        )

        fig.show()

        return fig
    return (CI_histograms,)


@app.cell
def _(CE_histograms):
    # Estimated effects
    def select_the_CE_of_choice(
        WPC_Relations, sel_display, sel_value, data, facility
    ):

        if sel_value in sel_display:

            position = WPC_Relations.index(sel_value)
            short_label = WPC_Relations[position]
            estimation_histogram = CE_histograms(data, sel_value, facility)

        else:
            estimation_histogram = None

        return estimation_histogram
    return (select_the_CE_of_choice,)


@app.function
def put_DAG_on(sel_value, sel_display, WPC, scm):

    if sel_value in sel_display:

        position = WPC.index(sel_value)

        label = scm[position]

        return label

    return None


@app.cell
def _():
    ###########################################################################################################################################
    return


@app.cell
def _(mo):
    # Display variables to choose
    variable_selector = mo.ui.dropdown(
        options=["Precipitation", "Specific Humidity", "Temperature", "Aerosols"],
        value="Temperature",
        label="Select Variable",
    )
    return (variable_selector,)


@app.cell
def _(
    data_M1_1_min,
    data_M1_30_min,
    data_S2_1_min,
    data_S2_30_min,
    data_S3_1_min,
    data_S3_30_min,
    selected_date,
):
    filtered_M1_1_min_data = select_date(data_M1_1_min, selected_date)
    filtered_S2_1_min_data = select_date(data_S2_1_min, selected_date)
    filtered_S3_1_min_data = select_date(data_S3_1_min, selected_date)

    filtered_M1_30_min_data = select_date(data_M1_30_min, selected_date)
    filtered_S2_30_min_data = select_date(data_S2_30_min, selected_date)
    filtered_S3_30_min_data = select_date(data_S3_30_min, selected_date)
    return (
        filtered_M1_1_min_data,
        filtered_S2_1_min_data,
        filtered_S3_1_min_data,
    )


@app.cell
def _(
    data_M1_30_min,
    data_S2_30_min,
    data_S3_30_min,
    filtered_M1_1_min_data,
    filtered_S2_1_min_data,
    filtered_S3_1_min_data,
    mo,
    plot_data,
    variable_selector,
):
    # Definición de las pestañas con el menú desplegable
    tab1 = mo.vstack(
        [
            variable_selector,
            mo.ui.table(data_M1_30_min),
            plot_data(filtered_M1_1_min_data, variable_selector.value),
        ]
    )

    tab2 = mo.vstack(
        [
            variable_selector,
            mo.ui.table(data_S2_30_min),
            plot_data(filtered_S2_1_min_data, variable_selector.value),
        ]
    )

    tab3 = mo.vstack(
        [
            variable_selector,
            mo.ui.table(data_S3_30_min),
            plot_data(filtered_S3_1_min_data, variable_selector.value),
        ]
    )

    # Creación de las pestañas
    tabs_30_min_data = mo.ui.tabs(
        {"Urban data": tab1, "Rural data": tab2, "Bay data": tab3}
    )
    return (tabs_30_min_data,)


@app.cell
def _(
    data_M1_1_min,
    data_S2_1_min,
    data_S3_1_min,
    mo,
    plot_data,
    variable_selector,
):
    tab1_1min = mo.vstack(
        [
            variable_selector,
            mo.ui.table(data_M1_1_min),
            plot_data(data_M1_1_min.head(1000), variable_selector.value),
        ]
    )

    tab2_1min = mo.vstack(
        [
            variable_selector,
            mo.ui.table(data_S2_1_min),
            plot_data(data_S2_1_min.head(1000), variable_selector.value),
        ]
    )

    tab3_1min = mo.vstack(
        [
            variable_selector,
            mo.ui.table(data_S3_1_min),
            plot_data(data_S3_1_min.head(1000), variable_selector.value),
        ]
    )

    tabs_1_min_data = mo.ui.tabs(
        {"Urban data": tab1_1min, "Rural data": tab2_1min, "Bay data": tab3_1min}
    )
    return (tabs_1_min_data,)


@app.cell
def _(mo, tabs_1_min_data, tabs_30_min_data):
    multi_tabs = mo.ui.tabs(
        {
            "30 minute interval data": tabs_30_min_data,
            "1 minute interval data": tabs_1_min_data,
        }
    )
    return (multi_tabs,)


@app.cell
def _(mo):
    select_causal_relation = mo.ui.radio(
        options=[
            "CO2-->Temperature",
            "Aerosols-->Temperature",
            "Specific Humidity-->Temperature",
            "Specific Humidity-->Aerosols",
        ],
        label="Select Causal Relation:",
    )
    return


@app.cell
def _(Causal_Relations, mo):
    select_causal_relation_CI_30_min = mo.ui.radio(
        options=Causal_Relations,
        label="Select Causal Relation:",
    )
    return (select_causal_relation_CI_30_min,)


@app.cell
def _(mo):
    select_causal_relation_CI_1_min = mo.ui.radio(
        options=[
            "Precipitation-->Precipitation",
            "Precipitation-->Temperature",
            "Precipitation-->Specific Humidity",
            "Precipitation-->Aerosols",
            "Temperature-->Precipitation",
            "Temperature-->Temperature",
            "Temperature-->Specific Humidity",
            "Temperature-->Aerosols",
            "Specific Humidity-->Precipitation",
            "Specific Humidity-->Temperature",
            "Specific Humidity-->Specific Humidity",
            "Specific Humidity-->Aerosols",
            "Aerosols-->Precipitation",
            "Aerosols-->Temperature",
            "Aerosols-->Specific Humidity",
            "Aerosols-->Aerosols",
        ],
        label="Select Causal Relation:",
    )
    return (select_causal_relation_CI_1_min,)


@app.cell
def _(
    M1_CD_coefs_1_min,
    M1_CI_coefs_30_min_m,
    M1_CI_coefs_30_min_n,
    S2_CD_coefs_1_min,
    S2_CI_coefs_30_min_m,
    S2_CI_coefs_30_min_n,
    S3_CD_coefs_1_min,
    S3_CI_coefs_30_min_m,
    S3_CI_coefs_30_min_n,
    select_causal_relation_CI_1_min,
    select_causal_relation_CI_30_min,
    select_the_CI_of_choice,
):
    M1_CI_coefs_30_min_stack_m = select_the_CI_of_choice(
        select_causal_relation_CI_30_min.value, M1_CI_coefs_30_min_m, "Urban", 1
    )
    S2_CI_coefs_30_min_stack_m = select_the_CI_of_choice(
        select_causal_relation_CI_30_min.value, S2_CI_coefs_30_min_m, "Rural", 1
    )
    S3_CI_coefs_30_min_stack_m = select_the_CI_of_choice(
        select_causal_relation_CI_30_min.value, S3_CI_coefs_30_min_m, "Bay", 1
    )

    M1_CI_coefs_30_min_stack_n = select_the_CI_of_choice(
        select_causal_relation_CI_30_min.value, M1_CI_coefs_30_min_n, "Urban", 1
    )
    S2_CI_coefs_30_min_stack_n = select_the_CI_of_choice(
        select_causal_relation_CI_30_min.value, S2_CI_coefs_30_min_n, "Rural", 1
    )
    S3_CI_coefs_30_min_stack_n = select_the_CI_of_choice(
        select_causal_relation_CI_30_min.value, S3_CI_coefs_30_min_n, "Bay", 1
    )

    M1_CD_coefs_1_min_stack = select_the_CI_of_choice(
        select_causal_relation_CI_1_min.value, M1_CD_coefs_1_min, "Urban", 5
    )

    S2_CD_coefs_1_min_stack = select_the_CI_of_choice(
        select_causal_relation_CI_1_min.value, S2_CD_coefs_1_min, "Rural", 5
    )

    S3_CD_coefs_1_min_stack = select_the_CI_of_choice(
        select_causal_relation_CI_1_min.value, S3_CD_coefs_1_min, "Bay", 5
    )
    return (
        M1_CD_coefs_1_min_stack,
        M1_CI_coefs_30_min_stack_m,
        M1_CI_coefs_30_min_stack_n,
        S2_CD_coefs_1_min_stack,
        S2_CI_coefs_30_min_stack_m,
        S2_CI_coefs_30_min_stack_n,
        S3_CD_coefs_1_min_stack,
        S3_CI_coefs_30_min_stack_m,
        S3_CI_coefs_30_min_stack_n,
    )


@app.cell
def _(
    CI_percentage_display,
    M1_CD_coefs_1_min_stack,
    M1_CI_percentage_1_min,
    S2_CD_coefs_1_min_stack,
    S2_CI_percentage_1_min,
    S3_CD_coefs_1_min_stack,
    S3_CI_percentage_1_min,
    mo,
    select_causal_relation_CI_1_min,
):
    carousel_CI_1_min_ = mo.carousel(
        [
            mo.hstack(
                [
                    M1_CD_coefs_1_min_stack,
                    CI_percentage_display(
                        select_causal_relation_CI_1_min.value,
                        M1_CI_percentage_1_min,
                        1,
                        10,
                    ),
                ]
            ),
            mo.hstack(
                [
                    S2_CD_coefs_1_min_stack,
                    CI_percentage_display(
                        select_causal_relation_CI_1_min.value,
                        S2_CI_percentage_1_min,
                        1,
                        10,
                    ),
                ]
            ),
            mo.hstack(
                [
                    S3_CD_coefs_1_min_stack,
                    CI_percentage_display(
                        select_causal_relation_CI_1_min.value,
                        S3_CI_percentage_1_min,
                        1,
                        10,
                    ),
                ]
            ),
        ]
    )

    rain_tabs = mo.ui.tabs(
        {
            "Precipitation Events": carousel_CI_1_min_,
        }
    )

    rain_CD_1_min_data = mo.vstack([select_causal_relation_CI_1_min, rain_tabs])
    return (rain_CD_1_min_data,)


@app.cell
def _(
    CI_percentage_display,
    M1_CI_coefs_30_min_stack_m,
    M1_CI_coefs_30_min_stack_n,
    M1_CI_percentage_30_min_m,
    S2_CI_coefs_30_min_stack_m,
    S2_CI_coefs_30_min_stack_n,
    S2_CI_percentage_30_min_m,
    S3_CI_coefs_30_min_stack_m,
    S3_CI_coefs_30_min_stack_n,
    S3_CI_percentage_30_min_m,
    mo,
    select_causal_relation_CI_30_min,
):
    carousel_CI_30_min_morning = mo.carousel(
        [
            mo.hstack(
                [
                    M1_CI_coefs_30_min_stack_m,
                    CI_percentage_display(
                        select_causal_relation_CI_30_min.value,
                        M1_CI_percentage_30_min_m,
                        1,
                        10,
                    ),
                ]
            ),
            mo.hstack(
                [
                    S2_CI_coefs_30_min_stack_m,
                    CI_percentage_display(
                        select_causal_relation_CI_30_min.value,
                        S2_CI_percentage_30_min_m,
                        1,
                        10,
                    ),
                ]
            ),
            mo.hstack(
                [
                    S3_CI_coefs_30_min_stack_m,
                    CI_percentage_display(
                        select_causal_relation_CI_30_min.value,
                        S3_CI_percentage_30_min_m,
                        1,
                        10,
                    ),
                ]
            ),
        ]
    )

    carousel_CI_30_min_night = mo.carousel(
        [
            mo.hstack(
                [
                    M1_CI_coefs_30_min_stack_n,
                    CI_percentage_display(
                        select_causal_relation_CI_30_min.value,
                        M1_CI_percentage_30_min_m,
                        1,
                        10,
                    ),
                ]
            ),
            mo.hstack(
                [
                    S2_CI_coefs_30_min_stack_n,
                    CI_percentage_display(
                        select_causal_relation_CI_30_min.value,
                        S2_CI_percentage_30_min_m,
                        1,
                        10,
                    ),
                ]
            ),
            mo.hstack(
                [
                    S3_CI_coefs_30_min_stack_n,
                    CI_percentage_display(
                        select_causal_relation_CI_30_min.value,
                        S3_CI_percentage_30_min_m,
                        1,
                        10,
                    ),
                ]
            ),
        ]
    )

    morning_night_tabs = mo.ui.tabs(
        {
            "Morning(4am-2pm)": carousel_CI_30_min_morning,
            "Night(4pm-2am)": carousel_CI_30_min_night,
        }
    )

    m_n_C1_30_min_data = mo.vstack(
        [select_causal_relation_CI_30_min, morning_night_tabs]
    )
    return (m_n_C1_30_min_data,)


@app.cell
def _(m_n_C1_30_min_data, mo, rain_CD_1_min_data):
    CI_30_min_tabs = mo.ui.tabs(
        {
            "1 minute interval data": rain_CD_1_min_data,
            "30 minute interval data": m_n_C1_30_min_data,
        }
    )
    return (CI_30_min_tabs,)


@app.cell
def _(WPC_Relations_1_min, mo):
    select_causal_relation_CE_1_min = mo.ui.radio(
        options=WPC_Relations_1_min, label="Select Causal Relation:"
    )
    return (select_causal_relation_CE_1_min,)


@app.cell
def _(WPC_Relations_30_min, mo):
    select_causal_relation_CE_30_min = mo.ui.radio(
        options=WPC_Relations_30_min, label="Select Causal Relation:"
    )
    return (select_causal_relation_CE_30_min,)


@app.cell
def _():
    ##################################### Causal Estimation ###############################
    return


@app.cell
def _(mo):
    scm_C02_Temp_diag = """
    %%{init: {'themeVariables': {'fontSize': '20px'}} }%%
    graph TB

    A((⠀⠀⠀⠀⠀CO2⠀⠀⠀⠀⠀))---> |30 minute lag|B((Temperature))
    C(("⠀Specific Humidity"))--> |30 minute lag|B((Temperature))
    D((Aerosols))--> |30 minute lag|B((⠀⠀Temperature⠀⠀))
    C--> |30 minute lag |D((⠀⠀⠀Aerosols⠀⠀⠀))

    style A fill:skyblue
    style B fill:skyblue
    style C fill:gray
    style D fill:gray


    linkStyle 0 stroke:red, stroke-width:3px
    linkStyle 1 stroke:gray, stroke-width:3px
    linkStyle 2 stroke:gray, stroke-width:3px
    linkStyle 3 stroke:gray, stroke-width:3px

    L_A_B_0@{ animation: fast}

    """
    md_scm_C02_Temp_diag = mo.mermaid(scm_C02_Temp_diag)

    ###################################################################################################################################################################

    scm_SH_Temp_diag = """
    %%{init: {'themeVariables': {'fontSize': '20px'}} }%%
    graph TB

    A((⠀⠀⠀⠀⠀CO2⠀⠀⠀⠀⠀))---> |30 minute lag|B((Temperature))
    C(("⠀Specific Humidity"))--> |30 minute lag|B((Temperature))
    D((Aerosols))--> |30 minute lag|B((⠀⠀Temperature⠀⠀))
    C--> |30 minute lag |D((⠀⠀⠀Aerosols⠀⠀⠀))

    style A fill:gray
    style B fill:skyblue
    style C fill:skyblue
    style D fill:gray


    linkStyle 0 stroke:gray, stroke-width:3px
    linkStyle 1 stroke:blue, stroke-width:3px
    linkStyle 2 stroke:gray, stroke-width:3px
    linkStyle 3 stroke:gray, stroke-width:3px

    L_C_B_0@{ animation: fast}

    """
    md_scm_SH_Temp_diag = mo.mermaid(scm_SH_Temp_diag)

    ###################################################################################################################################################################

    scm_aero_Temp_diag = """
    %%{init: {'themeVariables': {'fontSize': '20px'}} }%%
    graph TB

    A((⠀⠀⠀⠀⠀CO2⠀⠀⠀⠀⠀))---> |30 minute lag|B((Temperature))
    C(("⠀Specific Humidity"))--> |30 minute lag|B((Temperature))
    D((Aerosols))--> |30 minute lag|B((⠀⠀Temperature⠀⠀))
    C--> |30 minute lag |D((⠀⠀⠀Aerosols⠀⠀⠀))

    style A fill:gray
    style B fill:skyblue
    style C fill:gray
    style D fill:skyblue


    linkStyle 0 stroke:gray, stroke-width:3px
    linkStyle 1 stroke:gray, stroke-width:3px
    linkStyle 2 stroke:blue, stroke-width:3px
    linkStyle 3 stroke:gray, stroke-width:3px

    L_D_B_0@{ animation: fast}

    """
    md_scm_aero_Temp_diag = mo.mermaid(scm_aero_Temp_diag)

    ###################################################################################################################################################################


    scm_sh_aero_diag = """
    %%{init: {'themeVariables': {'fontSize': '20px'}} }%%
    graph TB

    A((⠀⠀⠀⠀⠀CO2⠀⠀⠀⠀⠀))---> |30 minute lag|B((Temperature))
    C(("⠀Specific Humidity"))--> |30 minute lag|B((Temperature))
    D((Aerosols))--> |30 minute lag|B((⠀⠀Temperature⠀⠀))
    C--> |30 minute lag |D((⠀⠀⠀Aerosols⠀⠀⠀))

    style A fill:gray
    style B fill:gray
    style C fill:skyblue
    style D fill:skyblue


    linkStyle 0 stroke:gray, stroke-width:3px
    linkStyle 1 stroke:gray, stroke-width:3px
    linkStyle 2 stroke:gray, stroke-width:3px
    linkStyle 3 stroke:blue, stroke-width:3px

    L_C_D_0@{ animation: fast}

    """
    md_scm_sh_aero_diag = mo.mermaid(scm_sh_aero_diag)

    ###################################################################################################################################################################

    scm_30_min = [
        md_scm_C02_Temp_diag,
        md_scm_aero_Temp_diag,
        md_scm_SH_Temp_diag,
        md_scm_sh_aero_diag,
    ]
    return scm_30_min, scm_sh_aero_diag


@app.cell
def _(mo):
    scm_SH_Temp_1_min_diag = """
    %%{init: {'themeVariables': {'fontSize': '20px'}} }%%

    graph LR
    A((⠀⠀Precipitation⠀⠀))
    B((⠀⠀Temperature⠀⠀))
    C((⠀Specific Humidity⠀))
    D((⠀⠀⠀Aerosols⠀⠀⠀⠀))


    C -->|5 minutes|B
    C -->|5 minutes|A
    C -->|5 minutes|D
    D --> |5 minutes|B
    A --> |instantaneous|D

    style A fill:gray
    style B fill:skyblue
    style C fill:skyblue
    style D fill:gray

    linkStyle 0 stroke:blue, stroke-width:3px
    linkStyle 1 stroke:gray, stroke-width:3px
    linkStyle 2 stroke:gray, stroke-width:3px
    linkStyle 3 stroke:gray, stroke-width:3px

    L_C_B_0@{ animation: fast}

    """
    md_scm_SH_Temp_1_min_diag = mo.mermaid(scm_SH_Temp_1_min_diag)

    ###################################################################################################################################################################

    scm_Precip_Aero_1_min_diag = """
    %%{init: {'themeVariables': {'fontSize': '20px'}} }%%

    graph LR
    A((⠀⠀Precipitation⠀⠀))
    B((⠀⠀Temperature⠀⠀))
    C((⠀Specific Humidity⠀))
    D((⠀⠀⠀Aerosols⠀⠀⠀⠀))


    C -->|5 minutes|B
    C -->|5 minutes|A
    C -->|5 minutes|D
    D --> |5 minutes|B
    A --> |instantaneous|D

    style A fill:skyblue
    style B fill:gray
    style C fill:gray
    style D fill:skyblue

    linkStyle 0 stroke:gray, stroke-width:3px
    linkStyle 1 stroke:gray, stroke-width:3px
    linkStyle 2 stroke:gray, stroke-width:3px
    linkStyle 3 stroke:gray, stroke-width:3px
    linkStyle 4 stroke:blue, stroke-width:3px

    L_A_D_0@{ animation: fast}

    """
    md_scm_Precip_Aero_1_min_diag = mo.mermaid(scm_Precip_Aero_1_min_diag)

    ###################################################################################################################################################################

    scm_SH_Aero_1_min_diag = """
    %%{init: {'themeVariables': {'fontSize': '20px'}} }%%

    graph LR
    A((⠀⠀Precipitation⠀⠀))
    B((⠀⠀Temperature⠀⠀))
    C((⠀Specific Humidity⠀))
    D((⠀⠀⠀Aerosols⠀⠀⠀⠀))


    C -->|5 minutes|B
    C -->|5 minutes|A
    C -->|5 minutes|D
    D --> |5 minutes|B
    A --> |instantaneous|D

    style A fill:gray
    style B fill:gray
    style C fill:skyblue
    style D fill:skyblue

    linkStyle 0 stroke:gray, stroke-width:3px
    linkStyle 1 stroke:gray, stroke-width:3px
    linkStyle 2 stroke:gray, stroke-width:3px
    linkStyle 3 stroke:blue, stroke-width:3px
    linkStyle 4 stroke:gray, stroke-width:3px

    L_C_D_0@{ animation: fast}

    """
    md_scm_SH_Aero_1_min_diag = mo.mermaid(scm_SH_Aero_1_min_diag)

    ###################################################################################################################################################################

    scm_SH_precip_1_min_diag = """
    %%{init: {'themeVariables': {'fontSize': '20px'}} }%%

    graph LR
    A((⠀⠀Precipitation⠀⠀))
    B((⠀⠀Temperature⠀⠀))
    C((⠀Specific Humidity⠀))
    D((⠀⠀⠀Aerosols⠀⠀⠀⠀))


    C -->|5 minutes|B
    C -->|5 minutes|A
    C -->|5 minutes|D
    D --> |5 minutes|B
    A --> |instantaneous|D

    style A fill:skyblue
    style B fill:gray
    style C fill:gray
    style D fill:skyblue

    linkStyle 0 stroke:gray, stroke-width:3px
    linkStyle 1 stroke:blue, stroke-width:3px
    linkStyle 2 stroke:gray, stroke-width:3px
    linkStyle 3 stroke:gray, stroke-width:3px
    linkStyle 4 stroke:gray, stroke-width:3px

    L_C_A_0@{ animation: fast}

    """
    md_scm_SH_precip_1_min_diag = mo.mermaid(scm_SH_precip_1_min_diag)

    ###################################################################################################################################################################

    scm_Aero_temp_1_min_diag = """
    %%{init: {'themeVariables': {'fontSize': '20px'}} }%%

    graph LR
    A((⠀⠀Precipitation⠀⠀))
    B((⠀⠀Temperature⠀⠀))
    C((⠀Specific Humidity⠀))
    D((⠀⠀⠀Aerosols⠀⠀⠀⠀))


    C -->|5 minutes|B
    C -->|5 minutes|A
    C -->|5 minutes|D
    D --> |5 minutes|B
    A --> |instantaneous|D

    style A fill:gray
    style B fill:skyblue
    style C fill:gray
    style D fill:skyblue

    linkStyle 0 stroke:gray, stroke-width:3px
    linkStyle 1 stroke:gray, stroke-width:3px
    linkStyle 2 stroke:gray, stroke-width:3px
    linkStyle 3 stroke:blue, stroke-width:3px
    linkStyle 4 stroke:gray, stroke-width:3px

    L_D_B_0@{ animation: fast}

    """
    md_scm_Aero_temp_1_min_diag = mo.mermaid(scm_Aero_temp_1_min_diag)

    ###################################################################################################################################################################

    scm_1_min = [
        md_scm_Precip_Aero_1_min_diag,
        md_scm_SH_Temp_1_min_diag,
        md_scm_SH_Aero_1_min_diag,
        md_scm_SH_precip_1_min_diag,
        md_scm_Aero_temp_1_min_diag,
    ]
    return (scm_1_min,)


@app.cell
def _(
    M1_WPC_1_min,
    M1_WPC_30_min_m,
    M1_WPC_30_min_n,
    S2_WPC_1_min,
    S2_WPC_30_min_m,
    S2_WPC_30_min_n,
    S3_WPC_1_min,
    S3_WPC_30_min_m,
    S3_WPC_30_min_n,
    WPC_Relations_1_min,
    WPC_Relations_30_min,
    select_causal_relation_CE_1_min,
    select_causal_relation_CE_30_min,
    select_the_CE_of_choice,
):
    M1_CE_1_min_stack = select_the_CE_of_choice(
        WPC_Relations_1_min,
        WPC_Relations_1_min,
        select_causal_relation_CE_1_min.value,
        M1_WPC_1_min,
        "Urban",
    )

    S2_CE_1_min_stack = select_the_CE_of_choice(
        WPC_Relations_1_min,
        WPC_Relations_1_min,
        select_causal_relation_CE_1_min.value,
        S2_WPC_1_min,
        "Rural",
    )

    S3_CE_1_min_stack = select_the_CE_of_choice(
        WPC_Relations_1_min,
        WPC_Relations_1_min,
        select_causal_relation_CE_1_min.value,
        S3_WPC_1_min,
        "Bay",
    )

    ################################################################

    M1_CE_30_min_stack_m = select_the_CE_of_choice(
        WPC_Relations_30_min,
        WPC_Relations_30_min,
        select_causal_relation_CE_30_min.value,
        M1_WPC_30_min_m,
        "Urban",
    )

    S2_CE_30_min_stack_m = select_the_CE_of_choice(
        WPC_Relations_30_min,
        WPC_Relations_30_min,
        select_causal_relation_CE_30_min.value,
        S2_WPC_30_min_m,
        "Rural",
    )

    S3_CE_30_min_stack_m = select_the_CE_of_choice(
        WPC_Relations_30_min,
        WPC_Relations_30_min,
        select_causal_relation_CE_30_min.value,
        S3_WPC_30_min_m,
        "Bay",
    )

    ################################################################

    M1_CE_30_min_stack_n = select_the_CE_of_choice(
        WPC_Relations_30_min,
        WPC_Relations_30_min,
        select_causal_relation_CE_30_min.value,
        M1_WPC_30_min_n,
        "Urban",
    )

    S2_CE_30_min_stack_n = select_the_CE_of_choice(
        WPC_Relations_30_min,
        WPC_Relations_30_min,
        select_causal_relation_CE_30_min.value,
        S2_WPC_30_min_n,
        "Rural",
    )

    S3_CE_30_min_stack_n = select_the_CE_of_choice(
        WPC_Relations_30_min,
        WPC_Relations_30_min,
        select_causal_relation_CE_30_min.value,
        S3_WPC_30_min_n,
        "Bay",
    )
    return (
        M1_CE_1_min_stack,
        M1_CE_30_min_stack_m,
        M1_CE_30_min_stack_n,
        S2_CE_1_min_stack,
        S2_CE_30_min_stack_m,
        S2_CE_30_min_stack_n,
        S3_CE_1_min_stack,
        S3_CE_30_min_stack_m,
        S3_CE_30_min_stack_n,
    )


@app.cell
def _(
    M1_CE_1_min_stack,
    S2_CE_1_min_stack,
    S3_CE_1_min_stack,
    WPC_Relations_1_min,
    mo,
    scm_1_min,
    select_causal_relation_CE_1_min,
):
    carousel_CE_1_min_ = mo.carousel(
        [
            mo.hstack(
                [
                    M1_CE_1_min_stack,
                    put_DAG_on(
                        select_causal_relation_CE_1_min.value,
                        WPC_Relations_1_min,
                        WPC_Relations_1_min,
                        scm_1_min,
                    ),
                ]
            ),
            mo.hstack(
                [
                    S2_CE_1_min_stack,
                    put_DAG_on(
                        select_causal_relation_CE_1_min.value,
                        WPC_Relations_1_min,
                        WPC_Relations_1_min,
                        scm_1_min,
                    ),
                ]
            ),
            mo.hstack(
                [
                    S3_CE_1_min_stack,
                    put_DAG_on(
                        select_causal_relation_CE_1_min.value,
                        WPC_Relations_1_min,
                        WPC_Relations_1_min,
                        scm_1_min,
                    ),
                ],
            ),
        ]
    )

    CE_rain_tabs = mo.ui.tabs(
        {
            "Precipitation Events": carousel_CE_1_min_,
        }
    )

    CE_rain_1_min_data = mo.vstack([select_causal_relation_CE_1_min, CE_rain_tabs])
    return (CE_rain_1_min_data,)


@app.cell
def _(
    M1_CE_30_min_stack_m,
    M1_CE_30_min_stack_n,
    S2_CE_30_min_stack_m,
    S2_CE_30_min_stack_n,
    S3_CE_30_min_stack_m,
    S3_CE_30_min_stack_n,
    WPC_Relations_30_min,
    mo,
    scm_30_min,
    select_causal_relation_CE_30_min,
):
    carousel_CE_30_min_m = mo.carousel(
        [
            mo.hstack(
                [
                    put_DAG_on(
                        select_causal_relation_CE_30_min.value,
                        WPC_Relations_30_min,
                        WPC_Relations_30_min,
                        scm_30_min,
                    ),
                    M1_CE_30_min_stack_m,
                ]
            ),
            mo.hstack(
                [
                    put_DAG_on(
                        select_causal_relation_CE_30_min.value,
                        WPC_Relations_30_min,
                        WPC_Relations_30_min,
                        scm_30_min,
                    ),
                    S2_CE_30_min_stack_m,
                ]
            ),
            mo.hstack(
                [
                    put_DAG_on(
                        select_causal_relation_CE_30_min.value,
                        WPC_Relations_30_min,
                        WPC_Relations_30_min,
                        scm_30_min,
                    ),
                    S3_CE_30_min_stack_m,
                ]
            ),
        ]
    )

    carousel_CE_30_min_n = mo.carousel(
        [
            mo.hstack(
                [
                    put_DAG_on(
                        select_causal_relation_CE_30_min.value,
                        WPC_Relations_30_min,
                        WPC_Relations_30_min,
                        scm_30_min,
                    ),
                    M1_CE_30_min_stack_n,
                ]
            ),
            mo.hstack(
                [
                    put_DAG_on(
                        select_causal_relation_CE_30_min.value,
                        WPC_Relations_30_min,
                        WPC_Relations_30_min,
                        scm_30_min,
                    ),
                    S2_CE_30_min_stack_n,
                ]
            ),
            mo.hstack(
                [
                    put_DAG_on(
                        select_causal_relation_CE_30_min.value,
                        WPC_Relations_30_min,
                        WPC_Relations_30_min,
                        scm_30_min,
                    ),
                    S3_CE_30_min_stack_n,
                ]
            ),
        ]
    )

    CE_morning_night_tabs = mo.ui.tabs(
        {
            "Morning(4am-2pm)": carousel_CE_30_min_m,
            "Night(4pm-2am)": carousel_CE_30_min_n,
        }
    )

    CE_m_n_30_min_data = mo.vstack(
        [select_causal_relation_CE_30_min, CE_morning_night_tabs]
    )
    return (CE_m_n_30_min_data,)


@app.cell
def _(CE_m_n_30_min_data, CE_rain_1_min_data, mo):
    CE_tabs = mo.ui.tabs(
        {
            "1 minute interval data": CE_rain_1_min_data,
            "30 minute interval data": CE_m_n_30_min_data,
        }
    )
    return (CE_tabs,)


@app.cell
def _(mo):
    Courage_display = mo.vstack(
        [
            mo.center(
                mo.md(
                    r"""
    <details>
    <summary>
        <span style="font-family:Times, monospace; font-size: 30px; color:deepskyblue">
            <center> Coastal-Urban-Rural Atmospheric Gradient Experiment</i> </center>
        </span>
    </summary>

    Launch: 
    <span style="font-size: 14px; color:blue"> 
    December 1st,2024 <br>
    </span>

    Principal Investigator: 
    <span style="font-size: 14px; color:blue">
    Kenneth Davis <br>
    </span>

    Team of 27 co-investigators <br>

    <center> 
      <hr>
      - Urban (M1)<br>
      Location:
      <span style="font-size: 14px; color:blue">
      Near Football field of the Heritage Early Learning Center building, Mayland, Baltimore.
       </span>

      <hr>
      Rural (S2)<br>
      Location: 
      <span style="font-size: 14px; color:blue">
      Frederick County, Maryland, about 35 miles northwest of Baltimore.
      </span>

      <hr>
      Bay (S3)<br>
      Location: 
      <span style="font-size: 14px; color:blue">
      Chesapeake Bay, Maryland, about 34 miles from Baltimore.
      </span>

    </center>

    </details>
    """
                )
            ),
            mo.center(
                mo.image(
                    src="Images/Baltimore_CoURAGE_Facility_IMAGE_jpeg.png",
                    width=600,
                    height=500,
                    rounded=True,
                )
            ),
        ]
    )
    return (Courage_display,)


@app.cell
def _():
    flowchart_diag = """
    %%{init: {'themeVariables': {'fontSize': '20px'}} }%%
    flowchart LR

        A@{ shape: cyl, label: "1 min. data" }
        B@{ shape: cyl, label: "30 min. data" }
        C["Meteorological variables"]
        D["Precipitation events"]
        E["Daily periods"]
        F[ <i>Tigramite </i>
        <p style="color:red;">Causal Discovery </p> to obtain correlation <br>coefficients using <br>PCMCI algorithm]
        G["Create SCM of physical interactions from the distribution of the cofficients "]
        H[<i>Tigramite</i> 
        <p style="color:red;">Causal Estimation </p>Wright's path analysis <br>doing interventions <br>to the causal variables]

        A --> C
        B --> C

        C --> |Temperature</>
               Specific Humidity</>
               Precipitation</>
               Aerosols| D

        C --> |Temperature</>
               Specific Humidity</>
               Precipitation</> 
               CO2</>
               Aerosols| E

        D --> |Event duration:</>
               >= 20 minutes</>
               Acceptable interuptions:</>
               <= 10 minutes| F

        E --> |Morning Period:</>
               4am - 2 pm</>
               Night Period:</>
               4pm - 2 am</>|F

       F-->G

       G-->H

        style A font-size:20px;
        style B font-size:20px;
        style C font-size:20px;
        style D font-size:20px;
        style E font-size:20px;
        style F font-size:20px;
        style G font-size:20px;
        style H font-size:20px;"""
    return (flowchart_diag,)


@app.cell
def _(flowchart_diag, mo):
    flowchart_diag_plot = mo.mermaid(flowchart_diag)
    return (flowchart_diag_plot,)


@app.cell
def _(Courage_display, date_plot, mo, multi_tabs, selected_date):
    C_data_display = mo.vstack(
        [Courage_display, date_plot, selected_date, multi_tabs]
    )
    return (C_data_display,)


@app.cell
def _(C_data_display, mo):
    mo.accordion(
        {
            "CoURAGE data": C_data_display,
        }
    )
    return


@app.cell
def _(flowchart_diag_plot, mo):
    mo.accordion(
        {
            "Project Methodology": flowchart_diag_plot,
        }
    )
    return


@app.cell
def _(CI_30_min_tabs, mo):
    mo.accordion({"Causal Inference": CI_30_min_tabs})
    return


@app.cell
def _(CI_30_min_tabs, mo):
    mo.accordion({"Causal Discovery": CI_30_min_tabs})
    return


@app.cell
def _(CE_tabs, mo):
    mo.accordion(
        {
            "Causal Estimation": CE_tabs,
        }
    )
    return


@app.cell
def _(mo, scm_sh_aero_diag):
    cm_diag_plot = mo.mermaid(scm_sh_aero_diag)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
