import pandas as pd

# Read the mature.csv and mirna_filtered.csv files
csv_file_mature = "mature.csv"
csv_file_filtered = "mirna_filtered.csv"

df_mature = pd.read_csv(csv_file_mature)
df_filtered = pd.read_csv(csv_file_filtered)

# Find the common data based on all columns
common_df = pd.merge(df_mature, df_filtered)

# Save the common data to a new CSV file
csv_file_common = "common.csv"
common_df.to_csv(csv_file_common, index=False)
print("Common data saved to common.csv")
