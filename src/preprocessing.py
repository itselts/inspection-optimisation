import pandas as pd
import geopandas as gpd
from shapely import wkt
from copy import deepcopy
import folium 
import os
import argparse
from points_curve import combined, plot_curve
import datetime

# Folium resources
# https://python-visualization.github.io/folium/quickstart.html
# https://python-visualization.github.io/folium/quickstart.html#GeoJSON/TopoJSON-Overlays 
# https://leafletjs.com/examples/geojson/ 

print(f"{'-'*40}")
print(f"{'Preprocessing...' : ^40}")


def build_cmd_args() -> argparse.Namespace:
    # Parser
    parser = argparse.ArgumentParser(description="Execute the whole work flow.")

    parser.add_argument('--request_ID', action='store', type=str, required=True)
    parser.add_argument('--s_date', action='store', required=True, help='Starting date of the horizon in iso format (yyyy-mm-dd)')
    parser.add_argument('--s_time', action='store', default='00:00', help='Optimisation start time in %H:%M format (e.g 07:00)')

    args = parser.parse_args()
    return args


args = build_cmd_args()
request_ID = args.request_ID
s_date_str = args.s_date
s_date = datetime.datetime(int(s_date_str[:4]), int(s_date_str[5:7]), int(s_date_str[8:10]), 7, 30)


# request_ID = "6E1C8D57-F535-4A71-AEC2-1D3AA4828DFA"
# s_date_str = "2023-06-16"
# s_date = datetime.datetime(int(s_date_str[:4]), int(s_date_str[5:7]), int(s_date_str[8:10]), 7, 30)


if not os.path.isdir(f"../runs/{request_ID}/outputs"):
    os.mkdir(f"../runs/{request_ID}/outputs")

##### DATA #####
# Inspections
insp_df = pd.read_parquet(f"../runs/{request_ID}/data/v_GetOptimiseInspectionInputItems.parquet")
insp_df['GeomWKT'] = insp_df['GeomWKT'].apply(wkt.loads) # Returns a geometric object from a wkt representation
insp_gdf = gpd.GeoDataFrame(insp_df, crs='epsg:4326', geometry='GeomWKT') # Geodataframe with WKT column as a MultiLineString
print(len(insp_df), "inspection segments.")

# Defects
def_df = pd.read_parquet(f"../runs/{request_ID}/data/v_GetOptimiseInputItems.parquet")
print(len(def_df), "defect jobs.")

# Depots
depot_df = pd.read_parquet(f"../runs/{request_ID}/data/v_DomainValueDepot.parquet")
print(len(depot_df), "depots.")

# Crews
crew_df =  pd.read_parquet(f"../runs/{request_ID}/data/v_GetOptimiseInputCrews.parquet")
print(len(crew_df), "crews.")


##### CREATING POINTS (Used in plot_solution.py) #####
# Defects
new_def_df = deepcopy(def_df.loc[:, ["ItemIdentifier", "LongStart", "LatStart"]])
new_def_df.rename({"ItemIdentifier": "ID", "LongStart": "X", "LatStart": "Y"}, axis=1, inplace=True)
new_def_df.ID = [ID+"0" for ID in new_def_df["ID"]]
new_def_df.reset_index(drop=True, inplace=True)
new_def_df.to_csv(f"../runs/{request_ID}/outputs/def_points.csv")

# Inspections
insp_list = []
for i, row in insp_df.iterrows():
    insp_list.append([''.join([str(row["ItemIdentifier"]), "0"]), row["LongStart"], row["LatStart"]])
    insp_list.append([''.join([str(row["ItemIdentifier"]), "1"]), row["LongEnd"], row["LatEnd"]])

new_insp_df = pd.DataFrame(insp_list)
if not insp_df.empty: # Check if it is a pure defect optimisation
    new_insp_df.columns = ["ID", "X", "Y"]
new_insp_df.to_csv(f"../runs/{request_ID}/outputs/insp_points.csv")

# Depots
new_depot_df = deepcopy(depot_df[["DepotRefNumber", "Long", "Lat"]])
new_depot_df.rename({"DepotRefNumber": "ID", "Long": "X", "Lat": "Y"}, axis=1, inplace=True)
new_depot_df.reset_index(drop=True, inplace=True)
new_depot_df.to_csv(f"../runs/{request_ID}/outputs/depot_points.csv")

##### COMBINING ALL POINTS #####
all_gdf = pd.concat([new_def_df, new_insp_df, new_depot_df])
all_gdf.reset_index(drop=True, inplace=True)
all_gdf.to_csv(f"../runs/{request_ID}/outputs/all_points.csv")



##### JOB SCORE #####
# Transforming the data (Setting the start day and time)
SCORES = deepcopy(def_df[["ItemIdentifier", "DueDate"]])
SCORES["OptStartDate"] = s_date
SCORES["DaysToDue"] = [(row["DueDate"] - s_date).total_seconds()/(60*60*24) for i, row in SCORES.iterrows()]

# Points curve for each job (Based on due date then scaled by activity multiplier)
SCORES['Score'] = SCORES.apply(lambda row: combined(
    x = -row['DaysToDue'], 
    growth = 1.25,
    due_score = 250,
    minimum_score = 50,
    cycle_height1 = 150, 
    cycle_maximum1 = 250,
    width1 = 10,
    cycle_height2 = 100, 
    cycle_maximum2 = 250,
    width2 = 5,
    act_multi = 1), axis=1)
SCORES.to_csv(f"../runs/{request_ID}/outputs/scores.csv", index=False)

plot_curve(request_ID) # Save score curve to outputs folder

##### PLOTTING #####
m = folium.Map()

if not insp_df.empty:
    x1,y1,x2,y2 = insp_gdf["GeomWKT"].total_bounds
    m.fit_bounds([[y1, x1], [y2, x2]])
else:
    x1,y1,x2,y2 = min(def_df['LongStart']), min(def_df['LatStart']), max(def_df['LongStart']), max(def_df['LatStart'])
    m.fit_bounds([[y1, x1], [y2, x2]])

# Inspections
for i, row in insp_gdf.iterrows():
    ASSET_ID = row["ItemIdentifier"]
    s, e = str(ASSET_ID) + "0", str(ASSET_ID) + "1"
    sx, sy = new_insp_df.loc[new_insp_df.ID == s, "X"].iloc[0], new_insp_df.loc[new_insp_df.ID == s, "Y"].iloc[0]
    ex, ey = new_insp_df.loc[new_insp_df.ID == e, "X"].iloc[0], new_insp_df.loc[new_insp_df.ID == e, "Y"].iloc[0]
    SEQ = row["DefaultSequence"]
    #DIR = row["Direction"]
    MLS = row["GeomWKT"]
    folium.GeoJson(MLS, tooltip = f"Asset ID: {ASSET_ID} \Default sequence: {SEQ}", style_function=lambda x: {"marker": "Circle", "color": "#0000FF", "weight": 5}).add_to(m)
    folium.CircleMarker(location=[float(sy), float(sx)], tooltip=f"ID = {s}", color="#FFA500", radius=3).add_to(m)
    folium.CircleMarker(location=[float(ey), float(ex)], tooltip=f"ID = {e}", color="#FFA500", radius=3).add_to(m)

# Defects
for i, row in def_df.iterrows():
    folium.CircleMarker(location=[row['LatStart'], row['LongStart']], tooltip=f'Defect ID: {row["ItemIdentifier"]}', color="#000000", radius=3).add_to(m)

# Depot
for i, row in depot_df.iterrows():
    folium.CircleMarker(location=[row['Lat'], row['Long']], tooltip=f'{row["dv_meaning"]} Depot', color="#FF0000", radius=5).add_to(m)

# Saving
m.save(f"../runs/{request_ID}/outputs/folium_map.html")

print("Geographical plot complete.")
print("Preprocessing complete.")
