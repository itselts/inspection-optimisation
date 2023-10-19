import pandas as pd
import geopandas as gpd
from shapely import wkt
import folium
import os
import numpy as np
from matplotlib.pyplot import cm
import argparse


# Folium resources
# https://python-visualization.github.io/folium/quickstart.html
# https://python-visualization.github.io/folium/quickstart.html#GeoJSON/TopoJSON-Overlays 
# https://leafletjs.com/examples/geojson/ 

print(f"{'-'*40}")
print(f"{'Creating solutions plot...' : ^40}")


def build_cmd_args() -> argparse.Namespace:
    # Parser
    parser = argparse.ArgumentParser(description="Execute the whole work flow.")

    parser.add_argument('--request_ID', action='store', type=str, required=True)
    #parser.add_argument('--s_date', action='store', required=True, help='Starting date of the horizon in iso format (yyyy-mm-dd)')
    #parser.add_argument('--s_time', action='store', default='00:00', help='Optimisation start time in %H:%M format (e.g 07:00)')

    args = parser.parse_args()
    return args


args = build_cmd_args()
request_ID = args.request_ID

##### DATA #####
# Inspections
insp_df = pd.read_parquet(f"../runs/{request_ID}/data/v_GetOptimiseInspectionInputItems.parquet")
insp_df['GeomWKT'] = insp_df['GeomWKT'].apply(wkt.loads) # Returns a geometric object from a wkt representation
insp_gdf = gpd.GeoDataFrame(insp_df, crs='epsg:4326', geometry='GeomWKT') # Geodataframe with WKT column as a MultiLineString

# Depots
depot_df = pd.read_parquet(f"../runs/{request_ID}/data/v_DomainValueDepot.parquet")

DEPOTS = list(depot_df["dv_code"])

INSPECTION_DATA = pd.read_csv(f"../runs/{request_ID}/outputs/insp_points.csv", header=0, index_col=0)
INSPECTIONS = [str(x) for x in list(INSPECTION_DATA.loc[:, "ID"])]

DEPOT_DATA = pd.read_csv(f"../runs/{request_ID}/outputs/depot_points.csv", header=0, index_col=0)
DEPOTS = [str(x) for x in list(DEPOT_DATA.loc[:, "ID"])]

NODES_DATA = pd.read_csv(f"../runs/{request_ID}/outputs/all_points.csv", header=0, index_col=0)
NODES_DATA.reset_index(drop=True, inplace=True)
NODES_DATA['ID'] = NODES_DATA['ID'].map(str)

##### RESULTS #####
i = 1
X_dict = {}
for f in os.listdir(f"../runs/{request_ID}/warm_start/X_matrices"):
    X = pd.read_csv(f"../runs/{request_ID}/warm_start/X_matrices/{f}", header=0)
    X.rename(dict(zip(range(len(X)), list(X.columns))), inplace=True)
    X_dict[i] = X
    i += 1

t = pd.read_csv(f"../runs/{request_ID}/warm_start/t.csv", index_col=False)
    

##### PLOTTING #####
##### CREW COLOURED PLOTS #####
m = folium.Map()
x1,y1,x2,y2 = insp_gdf["GeomWKT"].total_bounds
m.fit_bounds([[y1, x1], [y2, x2]])

n = len(X_dict)

colours = [tuple(rgb)[:3] for rgb in cm.rainbow(np.linspace(0, 1, n))]
colours = ['#{:02x}{:02x}{:02x}'.format(round(rgb[0]*255), round(rgb[1]*255), round(rgb[2]*255)) for rgb in colours] # In hex colours


for i, row in depot_df.iterrows():
    folium.CircleMarker(location=[row["Lat"], row["Long"]], color="#FF0000", radius=7).add_to(m)  # depot


for crew, X in X_dict.items():
    crew_color = colours[crew-1]
    for i, row in X.iterrows():
        for j, value in row.items():
            if round(value) == 1:

                # Plotting    
                if (i in INSPECTIONS) and (j in INSPECTIONS): # If both are segment start/end
                    if i[:-1] == j[:-1]: # Travelling on a segment (Blue)
                        segment = str(i[:-1])
                        mask = insp_gdf["ItemIdentifier"] == int(segment)
                        MLS = insp_gdf.loc[mask, "GeomWKT"].values[0]
                        REVERSIBLE = insp_gdf.loc[mask, "Reversible"].values[0]
                        folium.GeoJson(MLS, tooltip=f"Time: {t.loc[t['Identifier'] == int(j), 'Time'].values[0]} ID: {segment} Reversible: {REVERSIBLE}", style_function=lambda x: {"color": "#0000FF", "weight": 3}).add_to(m)
                        continue
                
                sx, sy = NODES_DATA.loc[NODES_DATA["ID"] == i, "X"].values[0], NODES_DATA.loc[NODES_DATA["ID"] == i, "Y"].values[0]
                ex, ey = NODES_DATA.loc[NODES_DATA["ID"] == j, "X"].values[0], NODES_DATA.loc[NODES_DATA["ID"] == j, "Y"].values[0]
                
                # Plotting the nodes
                if i in INSPECTIONS:
                    folium.CircleMarker(location=[float(sy), float(sx)], tooltip=f"ID = {i} coords: {sy, sx}", color="#FFA500", radius=3).add_to(m) 
                if j in INSPECTIONS:
                    folium.CircleMarker(location=[float(ey), float(ex)], tooltip=f"ID = {j} coords: {ey, ex}", color="#FFA500", radius=3).add_to(m)
                
                # Plotting the arcs
                locations = [(sy, sx), (ey, ex)]
                
                if (j not in DEPOTS): # Travelling between inspections or jobs
                    folium.PolyLine(locations=locations, opacity=0.6, color=crew_color, tooltip=f"Start: {i} End: {j} Time: {t.loc[t['Identifier'] == int(j), 'Time'].values[0]}").add_to(m)
                else: # Return to depot arc 
                    folium.PolyLine(locations=locations, opacity=0.6, color=crew_color).add_to(m)

                



# Saving 
m.save(f"../runs/{request_ID}/warm_start/folium_ws.html")

print("Solutions plot complete.")
