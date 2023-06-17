import pandas as pd

# Read the CSV file into a pandas DataFrame
df = pd.read_csv('common.csv')

# Apply the first filter to select desired columns
filtered_df = df[['Name', 'Biomarker_Name', 'Sequence', 'Sampling_Method', 'Collection_site', 'Disease_name', 'Regulation', 'Exosomal/Non-exosomal']]

# Apply the second filter to keep rows with the same 'Name' and 'Biomarker_Name'
filtered_df = filtered_df[filtered_df['Name'] == filtered_df['Biomarker_Name']]

# Apply the third filter to remove the 'Name' column
filtered_df = filtered_df.drop('Name', axis=1)

# Keep only one row for each unique 'Biomarker_Name'
filtered_df = filtered_df.drop_duplicates(subset='Biomarker_Name')

# Save the final filtered DataFrame to a new CSV file
filtered_df.to_csv('final_filtered_common.csv', index=False)
