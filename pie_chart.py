import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Load data from CSV files for different taxonomy levels
Phylum_data = pd.read_csv("Phylum-level.csv", delimiter=",")
Genus_data = pd.read_csv("Genus-level.csv", delimiter=",")
Species_data = pd.read_csv("Species-level.csv", delimiter=",")

# Calculate the overall abundance for each taxonomic level
Phylum_data['Overall_Abundance'] = Phylum_data.iloc[:, 1:].sum(axis=1)
Genus_data['Overall_Abundance'] = Genus_data.iloc[:, 1:].sum(axis=1)
Species_data['Overall_Abundance'] = Species_data.iloc[:, 1:].sum(axis=1)

# Sort the data by overall abundance
Phylum_data_sorted = Phylum_data.sort_values(by='Overall_Abundance', ascending=False)
Genus_data_sorted = Genus_data.sort_values(by='Overall_Abundance', ascending=False)
Species_data_sorted = Species_data.sort_values(by='Overall_Abundance', ascending=False)

# Custom label formatting function to show only percentages on pie chart
def custom_autopct(pct, allvals):
    """Custom label formatting for pie chart slices."""
    return f"{pct:.1f}%" if pct >= 1 else ''

# Function to create a pie chart showing the abundance of taxa
def pie_chart_showing_abundance(df, figsize=(10, 7), dpi=300, cmap_name='gist_earth', Taxa="Taxa"):
    """Generates and displays a pie chart showing the overall abundance of taxa."""
    # Create a figure object
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)
    
    # Generate color mapping from the specified colormap
    cmap = plt.get_cmap(cmap_name)
    colors = cmap(np.linspace(0, 1, len(df["Overall_Abundance"])))
    
    # Creating the pie chart using the "Overall_Abundance" column directly
    wedges, texts, autotexts = ax.pie(
        df["Overall_Abundance"],
        startangle=90,
        colors=colors,  # Use colors from the colormap
        autopct=lambda pct: custom_autopct(pct, df["Overall_Abundance"]),
        pctdistance=0.75  # Adjust this to move text closer or farther from the center
    )

    # Adjust the position of each text label to match the angle of the wedge
    for text, autotext in zip(texts, autotexts):
        ang = (autotext.get_position()[0]**2 + autotext.get_position()[1]**2)**0.5
        x = autotext.get_position()[0]
        y = autotext.get_position()[1]
        if x < 0:
            ang = np.degrees(np.arctan2(y, x))
            autotext.set_rotation(ang - 180)
            autotext.set_color('white')
        else:
            ang = np.degrees(np.arctan2(y, x))
            autotext.set_rotation(ang)
            autotext.set_color('white')
    
    # Ensure the pie chart is drawn as a circle
    ax.axis('equal')
    
    # Set the title of the chart
    plt.title(f"Overall Abundance of {Taxa} in Group A")
    
    # Create a legend with taxonomic names, using the colors of the slices
    plt.legend(
        wedges,
        df[Taxa],
        title=f"{Taxa} Names",
        loc="center right",
        bbox_to_anchor=(0.9, 0, 0.5, 1),
        fontsize=7,
        ncol=2
    )
    
    # Adjust layout to not cut off content
    plt.tight_layout()
    
    # Display the plot
    plt.show()

# Generate the pie chart for each taxonomic level
pie_chart_showing_abundance(Phylum_data_sorted, Taxa="Phylum")
pie_chart_showing_abundance(Genus_data_sorted, Taxa="Genus")
pie_chart_showing_abundance(Species_data_sorted, Taxa="Species")
