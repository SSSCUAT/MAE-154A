import pandas as pd
import matplotlib.pyplot as plt
import os

def parse_apc_file(filepath):
    """
    Parses APC propeller performance text files into a Pandas DataFrame.
    """
    data = []
    current_rpm = None
    read_data = False
    
    # Extract the propeller name from the filename (e.g., 'PER3_16x7')
    prop_name = os.path.basename(filepath).replace('.txt', '')
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Identify the RPM block start
            if "PROP RPM =" in line:
                try:
                    current_rpm = float(line.split("=")[1].strip())
                    read_data = False # Reset until we hit the table headers
                except ValueError:
                    pass
                continue
                
            # The actual numerical data starts right after the unit headers
            if "(mph)" in line and "(Adv_Ratio)" in line:
                read_data = True
                continue
                
            # A blank line indicates the end of the current RPM's data block
            if read_data and not line:
                read_data = False
                continue
                
            # Parse the numeric data lines
            if read_data:
                parts = line.split()
                # Ensure we are reading a valid data line (APC files usually have 15 columns)
                if len(parts) >= 15:
                    try:
                        row = {
                            'Propeller': prop_name,
                            'RPM': current_rpm,
                            'Speed V (mph)': float(parts[0]),
                            'Advance Ratio (J)': float(parts[1]),
                            'Efficiency (Pe)': float(parts[2]),
                            'Thrust (Lbf)': float(parts[7]),
                            'Power (W)': float(parts[8])
                        }
                        data.append(row)
                    except ValueError:
                        pass # Skips any unexpected text lines

    return pd.DataFrame(data)

# 1. Load the data
# Make sure the files 'PER3_16x7.txt' and 'PER3_16x8.txt' are in the same folder as this script
files_to_process = ['PER3_16x7.txt', 'PER3_16x8.txt']
dataframes = []

for file in files_to_process:
    if os.path.exists(file):
        dataframes.append(parse_apc_file(file))
    else:
        print(f"File not found: {file}")

if not dataframes:
    raise ValueError("No data files were found to process.")

# Combine into a single DataFrame
df_all = pd.concat(dataframes, ignore_index=True)

# 2. Select an RPM to analyze 
# Efficiency curves change drastically with RPM. 
# We'll plot at 5000 RPM (a typical cruising RPM for 16-inch RC props).
target_rpm = 1000.0
df_rpm = df_all[df_all['RPM'] == target_rpm]

if df_rpm.empty:
    print(f"No data found for {target_rpm} RPM. Available RPMs: {df_all['RPM'].unique()}")
else:
    # 3. Create Plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'Propeller Performance Comparison at {int(target_rpm)} RPM', fontsize=16)

    # Plot 1: Efficiency vs. Advance Ratio (J)
    # The peak of this curve tells you the aerodynamic limits of the prop.
    for prop in df_rpm['Propeller'].unique():
        subset = df_rpm[df_rpm['Propeller'] == prop]
        axes[0].plot(subset['Advance Ratio (J)'], subset['Efficiency (Pe)'], marker='.', label=prop)
    
    axes[0].set_title('Efficiency vs. Advance Ratio')
    axes[0].set_xlabel('Advance Ratio (J)')
    axes[0].set_ylabel('Efficiency (Pe)')
    axes[0].grid(True)
    axes[0].legend()

    # Plot 2: Efficiency vs. Speed (mph)
    # This helps you match the propeller to your RC plane's expected cruise speed.
    for prop in df_rpm['Propeller'].unique():
        subset = df_rpm[df_rpm['Propeller'] == prop]
        axes[1].plot(subset['Speed V (mph)'], subset['Efficiency (Pe)'], marker='.', label=prop)
    
    axes[1].set_title('Efficiency vs. Airspeed')
    axes[1].set_xlabel('Speed (mph)')
    axes[1].set_ylabel('Efficiency (Pe)')
    axes[1].grid(True)
    axes[1].legend()

    # Plot 3: Thrust vs. Speed (mph)
    # Important to ensure the prop produces enough thrust at your target speed.
    for prop in df_rpm['Propeller'].unique():
        subset = df_rpm[df_rpm['Propeller'] == prop]
        axes[2].plot(subset['Speed V (mph)'], subset['Thrust (Lbf)'], marker='.', label=prop)
    
    axes[2].set_title('Thrust vs. Airspeed')
    axes[2].set_xlabel('Speed (mph)')
    axes[2].set_ylabel('Thrust (Lbf)')
    axes[2].grid(True)
    axes[2].legend()

    plt.tight_layout()
    plt.show()