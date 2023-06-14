import pandas as pd

# Read the mirna.csv file
csv_file_mirna = "mirna.csv"
df_mirna = pd.read_csv(csv_file_mirna)

# Filter rows based on Disease_Category
df_filtered = df_mirna[df_mirna['Disease_Category'] != 'Healthy']

# Save the filtered DataFrame to a new CSV file
csv_file_filtered = "mirna_filtered.csv"
df_filtered.to_csv(csv_file_filtered, index=False)
print("Filtered data saved to mirna_filtered.csv")
