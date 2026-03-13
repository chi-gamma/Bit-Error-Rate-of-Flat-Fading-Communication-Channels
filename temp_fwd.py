import pandas as pd
import yfinance as yf
import numpy as np


def longest_monotonic_subseq_idx(values, decreasing=True):
    values = np.asarray(values, dtype=float)
    n = len(values)
    
    dp = np.ones(n, dtype=int)
    prev = np.full(n, -1)
    
    for i in range(n):
        for j in range(i):
            cond = values[j] >= values[i] if decreasing else values[j] <= values[i]
            if cond and dp[j] + 1 > dp[i]:
                dp[i] = dp[j] + 1
                prev[i] = j
                
    # reconstruct sequence
    idx = np.argmax(dp)
    seq = []
    while idx != -1:
        seq.append(idx)
        idx = prev[idx]
        
    return seq[::-1]


def clean_monotonic(df):
    df = df.sort_values(by="strike", ignore_index=True)
    
    flag = True if df["option_type"].iloc[0] == "call" else False
    keep = longest_monotonic_subseq_idx(df["mid"], decreasing=flag)
    
    return df.iloc[keep]


def remove_spread_outliers(df):
    spread = df['ask'] - df['bid']
    mu = spread.mean()
    sd = spread.std()
    if sd > 0:
        df = df[spread < mu + 3*sd]
    return df



def filter_data(df):
    """Clean and filter raw option chain data."""
    
    df = df.copy()
    
    # Parse dates
    date_cols = ["as_of_date", "expiry", "last_trade_date"]
    for col in date_cols:
        df[col] = pd.to_datetime(df[col], errors="coerce")
    df = df.dropna(subset=date_cols)
        
    # Drop options with very stale last trade date
    df = df[df["last_trade_date"] >= (df["as_of_date"] - pd.offsets.BusinessDay(n=3))] 
    
    # Require valid two-sided quotes
    df = df[(df["bid"] >= 0.001) & (df["ask"] >= 0.001)]
    
    df["mid"] = df[["bid", "ask"]].mean(axis=1)
    
    # Relative spread filter
    df = df[df["ask"] <= 3.0 * df["bid"]]
    
    df = df.groupby("expiry", group_keys=False).apply(remove_spread_outliers)
    df = df.groupby(["expiry", "option_type"], group_keys=False).apply(clean_monotonic)

    # Tenor
    df['tenor'] = np.round((df['expiry'] - df['as_of_date']).dt.days / 365.0, 6)
    
    return df.sort_values(['tenor', 'strike'], ignore_index=True)


# read discount curve
rates = pd.read_excel('C:/Users/nwozo/Downloads/Forwards_DFs.xlsx')
log_discounts = np.log(rates['DF'])
curve_times = rates['T'].to_numpy()
discount_curve = lambda t: np.exp(np.interp(t, curve_times, log_discounts))


# read and process option chain
options = pd.read_csv("C:/Users/nwozo/Documents/option_chain.csv")
options = filter_data(options)


grouped = options.groupby("expiry") 
for (expiry,), group in grouped:
    group = group.copy()
    
    T = group["tenor"].iloc[0]
    df = discount_curve(T)
    
    ["strike", "bid", "ask", "mid"]
    
    
    calls = group[(group["option_type"]=="call") & (group["expiry"]==expiry)][["strike", "bid", "ask", "mid"]].copy()
    puts = group[(group["option_type"]=="put") & (group["expiry"]==expiry)][["strike", "bid", "ask", "mid"]].copy()
    calls.columns = ["strike"] + ["call_" + col for col in calls.columns[1:]]
    puts.columns = ["strike"] + ["put_" + col for col in puts.columns[1:]]

    merged = pd.merge(calls, puts, on='strike')
    if len(merged) < 4:
        continue
    
    f_min = (merged["strike"] + (merged["call_bid"] - merged["put_ask"]) / df).to_numpy()
    f_max = (merged["strike"] + (merged["call_ask"] - merged["put_bid"]) / df).to_numpy()
    
    f_low, f_high = f_min.max(), f_max.min() # admissible range of forwards
    
    if f_low > f_high:
        continue
    
    merged['mid_fwd'] = merged["strike"] + (merged["call_mid"] - merged["put_mid"]) / df
    
    max_iter = 10
    j = 0
    f_guess = (f_low + f_high) / 2.0
    while j <= max_iter:
        data1 = merged.loc[[ merged[merged['strike'] < f_guess]['strike'].idxmax()]]
        data2 = merged.loc[[ merged[merged['strike'] > f_guess]['strike'].idxmin()]]
        combined = pd.concat([data1, data2], ignore_index=True)
        w = 0.0
        if len(combined)==0:
            continue
        elif len(combined)==1:
            new_f_guess = combined['mid_fwd'].iloc[0]
            f_guess = 0.3 * new_f_guess + 0.7 * f_guess
            continue
        else:
            new_f_guess = combined['mid_fwd'].mean()
            change = abs(new_f_guess - f_guess) / f_guess
            f_guess = new_f_guess
            k1, k2 = tuple(combined['strike'])
            if f_guess >= k1 and f_guess <= k2 and change <= 0.001:
                break
            
            
            
            

    
    
    def get_fwd(f_low, f_high, merged):
        f_guess = (f_low + f_high) / 2.0
        
        data = pd.concat([
                    merged.loc[[ merged[merged['strike'] > f_guess]['strike'].idxmin()]],
                    merged.loc[[ merged[merged['strike'] < f_guess]['strike'].idxmax()]]
                ], ignore_index=True)
        if len(data) != 2:
            continue
        
        new_f_guess = np.mean(data["strike"] + (data["call_mid"] - data["put_mid"]) / df)
    
    

    

    
    
max_iter = 10
j = 0
fwd = data['strike'].max() #np.mean(data['mid_f'])
while j <= max_iter:
    data1 = data[(data['strike'] < fwd) & (data['low_f'] < fwd) &
                 (data['high_f'] > fwd)].tail(1)
    data2 = data[(data['strike'] > fwd) & (data['low_f'] < fwd) &
                 (data['high_f'] > fwd)].head(1)
    if len(data1)==0 and len(data2)==0: break
    else:
        new_forward = pd.concat([data1, data2])['mid_f'].mean()
        change = abs(new_forward - fwd) / fwd
        fwd = new_forward
        if len(data1)==1 and len(data2)==1 and \
        (data2.iloc[0,0] < new_forward or data1.iloc[0,0] > new_forward) : continue
        if change < 0.001: break
    j += 1

df_fwd.loc[i, 'date'] = date
df_fwd.loc[i, 'tenor'] = T
df_fwd.loc[i, 'forward'] = fwd
df_fwd.loc[i, 'DF'] = D
i += 1



group = options[options['expiry']==options.loc[0, 'expiry']].copy()



abe = pd.concat([abc, abd])












        
        
        
        
        
        
        























obj = yf.Ticker("^SPX")
dates = obj.options

options = pd.DataFrame()
for date in dates:
    calls, puts, _ = obj.option_chain(date)
    calls['option_type'] = 'call'
    puts['option_type'] = 'put'
    df = pd.concat([calls, puts])
    df['expiry'] = date
    options = pd.concat([options, df])
    
options['last_trade_date'] = options['lastTradeDate'].dt.strftime("%Y-%m-%d")
options['as_of_date'] = '2026-03-12'
options = options[['as_of_date', 'expiry', 'last_trade_date', 'option_type', 'strike', 'bid', 'ask']]


options.to_csv("C:/Users/nwozo/Documents/option_chain_20260312.csv")









    
    
    
options['as_of_date'] = '2026-03-06'
options.columns
options['contractSize'].unique()
options = options[['option_type', 'as_of_date', 'exp_date', 'lastTradeDate', 'strike', 'bid', 'ask', 'volume']]
options.columns = ['option_type', 'as_of_date', 'exp_date', 'last_trade_date', 'strike', 'bid', 'ask', 'volume']
options.to_csv('C:/Users/nwozo/Documents/option_chain_20260306.csv')
options['last_trade_date'].dt.date.strftime('%Y-%m-%d')
options['last_trade_date'].dt.date
options['last_trade_date'].dt.date.dt.strftime('%Y-%m-%d')
abc = options['last_trade_date'].dt.date
abc.dt
abc[0]
options.index = np.arange(len(options))
import numpy as np
options.index = np.arange(len(options))
abc = options['last_trade_date'].dt.date
abc[0]
abc = pd.to_datettime(options['last_trade_date'].dt.date)
abc = pd.to_datetime(options['last_trade_date'].dt.date)
abc = pd.to_datetime(options['last_trade_date'].dt.date).dt.strftime('%Y-%m-%d')
options['last_trade_date'] = pd.to_datetime(options['last_trade_date'].dt.date).dt.strftime('%Y-%m-%d')
options.to_csv('C:/Users/nwozo/Documents/option_chain_20260306.csv')
%clear
import pandas as pd
import yfinance as yf
import numpy as np


def filter_data(df):
    """Clean and filter raw option chain data."""
    
    df = df.copy()
    
    # Require valid two-sided quotes
    df = df[(df['bid'] >= 0.05) & (df['ask'] >= 0.05)]
    
    # Relative spread filter
    rel_spread = (df['ask'] - df['bid']) / df['mid']
    df = df[rel_spread <= 1.0]
    
    # Parse dates
    date_cols = ['as_of_date', 'exp_date', 'last_trade_date']
    for col in date_cols:
        df[col] = pd.to_datetime(df[col], errors='coerce')
    
    df = df.dropna(subset=date_cols)
    
    # Drop options with very stale last trade date
    df = df[df['last_trade_date'] >= (df['as_of_date'] - pd.offsets.BusinessDay(n=3))] 
    
    # Tenor
    df['tenor'] = np.round((df['exp_date'] - df['as_of_date']).dt.days / 365.0, 6)
    
    # Remove spread outliers per tenor
    def remove_spread_outliers(group):
        spread = group['ask'] - group['bid']
        mu = spread.mean()
        sd = spread.std()
        if sd > 0:
            group = group[spread < mu + 3*sd]
        return group
    
    df = df.groupby('tenor', group_keys=False).apply(remove_spread_outliers)
    
    return df.sort_values(['tenor', 'strike'], ignore_index=True)


df = pd.read_csv('C:/Users/nwozo/Documents/option_chain_20260306.csv')
df2 = filter_data(df)
data = pd.read_csv('C:/Users/nwozo/Documents/option_chain_20260306.csv')
df = data.copy()

df['mid'] = 0.5 * (df['ask'] + df['bid'])

def filter_data(data):
    """Clean and filter raw option chain data."""
    
    df = data.copy()
    
    df['mid'] = 0.5 * (df['ask'] + df['bid'])
    
    # Require valid two-sided quotes
    df = df[(df['bid'] >= 0.05) & (df['ask'] >= 0.05)]
    
    # Relative spread filter
    rel_spread = (df['ask'] - df['bid']) / df['mid']
    df = df[rel_spread <= 1.0]
    
    # Parse dates
    date_cols = ['as_of_date', 'exp_date', 'last_trade_date']
    for col in date_cols:
        df[col] = pd.to_datetime(df[col], errors='coerce')
    
    df = df.dropna(subset=date_cols)
    
    # Drop options with very stale last trade date
    df = df[df['last_trade_date'] >= (df['as_of_date'] - pd.offsets.BusinessDay(n=3))] 
    
    # Tenor
    df['tenor'] = np.round((df['exp_date'] - df['as_of_date']).dt.days / 365.0, 6)
    
    # Remove spread outliers per tenor
    def remove_spread_outliers(group):
        spread = group['ask'] - group['bid']
        mu = spread.mean()
        sd = spread.std()
        if sd > 0:
            group = group[spread < mu + 3*sd]
        return group
    
    df = df.groupby('tenor', group_keys=False).apply(remove_spread_outliers)
    
    return df.sort_values(['tenor', 'strike'], ignore_index=True)
df = filter_data(data)
data = pd.read_excel('C:/Users/nwozo/Downloads/Forwards_DFs.xlsx')
ln_dfs = np.log(data['DF'])
df_interp = lambda t: np.exp(np.interp(t, data['T'], ln_dfs))


spot = 6740.02
res = pd.DataFrame()
dates = df['exp_date'].unique()
reg_cols = ['mid_call', 'mid_put', 'strike']
for i, date in enumerate(dates):
    T = df.loc[df['exp_date']==date, 'tenor'].values[0]
    calls = df[(df['option_type']=='call') & (df['exp_date']==date)].copy()
    puts = df[(df['option_type']=='put') & (df['exp_date']==date)].copy()
    calls['mid_call'] = calls['mid']
    puts['mid_put'] = puts['mid']
    merged = pd.merge(calls, puts, how='left', on='strike')
    merged = merged[reg_cols]
    merged = merged.dropna()
    merged['mid'] = merged['mid_call'] - merged['mid_put']
    if len(merged) < 6: continue
    merged = merged.sort_values(by='mid', key=abs, ignore_index=True) #[:4]
    df = df_interp(T)
    N = len(merged)
    X = spot * np.ones(N).reshape(-1, 1)
    y = (merged['strike'] + merged['mid'] / df).to_numpy() 
    beta = np.linalg.lstsq(X, y, rcond=None)[0]  
    fwd = spot * beta[0]
    res.loc[i, 'Date'] = date
    res.loc[i, 'T'] = T
    res.loc[i, 'DF'] = df
    res.loc[i, 'Fwd'] = fwd
df = filter_data(data)
data = pd.read_csv('C:/Users/nwozo/Documents/option_chain_20260306.csv')
df = filter_data(data)

data2 = pd.read_excel('C:/Users/nwozo/Downloads/Forwards_DFs.xlsx')
ln_dfs = np.log(data2['DF'])
df_interp = lambda t: np.exp(np.interp(t, data2['T'], ln_dfs))

spot = 6740.02
res = pd.DataFrame()
dates = df['exp_date'].unique()
reg_cols = ['mid_call', 'mid_put', 'strike']
for i, date in enumerate(dates):
    T = df.loc[df['exp_date']==date, 'tenor'].values[0]
    calls = df[(df['option_type']=='call') & (df['exp_date']==date)].copy()
    puts = df[(df['option_type']=='put') & (df['exp_date']==date)].copy()
    calls['mid_call'] = calls['mid']
    puts['mid_put'] = puts['mid']
    merged = pd.merge(calls, puts, how='left', on='strike')
    merged = merged[reg_cols]
    merged = merged.dropna()
    merged['mid'] = merged['mid_call'] - merged['mid_put']
    if len(merged) < 6: continue
    merged = merged.sort_values(by='mid', key=abs, ignore_index=True) #[:4]
    discount = df_interp(T)
    N = len(merged)
    X = spot * np.ones(N).reshape(-1, 1)
    y = (merged['strike'] + merged['mid'] / discount).to_numpy() 
    beta = np.linalg.lstsq(X, y, rcond=None)[0]  
    fwd = spot * beta[0]
    res.loc[i, 'Date'] = date
    res.loc[i, 'T'] = T
    res.loc[i, 'DF'] = discount
    res.loc[i, 'Fwd'] = fwd
q = -np.log((res['Fwd'] * res['DF']) / spot) / res['T']
res.loc[4, 'Fwd'] = 6742.031 
q = -np.log((res['Fwd'] * res['DF']) / spot) / res['T']
df['tenor'].round(8).unique()








def EstimateForwards(df, df_interp):
    options = df.copy()
    df_fwd = pd.DataFrame()
    groups = options.groupby('exp_date')
    i = 0
    for (date, group,) in groups:
        T = group['tenor'].values[0]
        D = df_interp(T) # interpolate discount factor
        calls = df[(options['exp_date']==date) & (options['option_type']=='call')].loc[:][['strike', 'bid', 'ask', 'mid']]
        puts = df[(options['exp_date']==date) & (options['option_type']=='put')].loc[:][['strike', 'bid', 'ask', 'mid']]
        data = pd.merge(calls, puts, on='strike').dropna()
        data['mid'] = data['mid_x'] - data['mid_y']
        data = data.sort_values(by='mid', key=abs, ignore_index=True)[:6]
        
        if len(data)==0: continue
        temp1 = data['strike'] + (data['bid_x'] - data['bid_y']) / D
        temp2 = data['strike'] + (data['ask_x'] - data['ask_y']) / D
        data['low_f'] = np.minimum(temp1, temp2)
        data['high_f'] = np.maximum(temp1, temp2)
        data['mid_f'] = data['strike'] + (data['mid_x'] - data['mid_y']) / D
        
        max_iter = 10
        j = 0
        fwd = data['strike'].max() #np.mean(data['mid_f'])
        while j <= max_iter:
            data1 = data[(data['strike'] < fwd) & (data['low_f'] < fwd) &
                         (data['high_f'] > fwd)].tail(1)
            data2 = data[(data['strike'] > fwd) & (data['low_f'] < fwd) &
                         (data['high_f'] > fwd)].head(1)
            if len(data1)==0 and len(data2)==0: break
            else:
                new_forward = pd.concat([data1, data2])['mid_f'].mean()
                change = abs(new_forward - fwd) / fwd
                fwd = new_forward
                if len(data1)==1 and len(data2)==1 and \
                (data2.iloc[0,0] < new_forward or data1.iloc[0,0] > new_forward) : continue
                if change < 0.001: break
            j += 1
        
        df_fwd.loc[i, 'date'] = date
        df_fwd.loc[i, 'tenor'] = T
        df_fwd.loc[i, 'forward'] = fwd
        df_fwd.loc[i, 'DF'] = D
        i += 1
    
    df_fwd = df_fwd.dropna()
    return df_fwd
