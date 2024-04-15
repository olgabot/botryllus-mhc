
def add_percent_rolling(df):
    reviewed_total = df["len"].sum()

    df["percent"] = 100 * df["len"] / reviewed_total
    df = df.sort_values("percent", ascending=False)
    df["percent_rolling"] = df.percent.cumsum()
    return df
