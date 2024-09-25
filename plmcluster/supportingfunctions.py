def apply_func_to_window(df, col, func, window_size):
    return df[col].rolling(window=window_size).apply(lambda x: func(x))