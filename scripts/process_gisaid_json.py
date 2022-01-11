# Load pandas library
import pandas as pd

# Metadata dataframe
metadata_df = pd.DataFrame()

# Select chunk size. Insert a lower value if process is killed by memory usage
chunk_size = 100000

# Load json in chunks
gisaid_metadata_chunk = pd.read_json("json_file.json", lines = True, chunksize = chunk_size)

# Json to dataframe splitting covv_location by " / " for each chunk 
count = 0
for chunk_df in gisaid_metadata_chunk:
	count += 1
	print("Loading chunk number " + str(count))

	max_slashes = max(chunk_df.covv_location.str.count(" / ")) # Get max number of slashes in covv_location field
	splited_location = chunk_df["covv_location"].str.split(" / ", n = max_slashes, expand = True) # Split covv_location by " / "
	for i in range(0, max_slashes + 1): # Add new columns
		column_name = "Division_" + str(i)
		chunk_df[column_name] = splited_location[i]
	chunk_df.drop(columns = ["covv_location"], inplace = True) # Remove covv_location column
	metadata_df = metadata_df.append(chunk_df)

# Just to check
usa = metadata_df[metadata_df["Division_1"] == "USA"]
sub_usa = usa.sample(n = 4000)

print(sub_usa)
print(metadata_df.columns)
print(metadata_df.index)
print(len(usa.index))
print(len(sub_usa.index))