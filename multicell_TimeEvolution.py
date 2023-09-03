import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from matplotlib.animation import FuncAnimation
from PIL import Image
import tempfile
import os

# Load the data from the text file
data = np.loadtxt("evol_temp_redA.txt")

# Extract time and concentration columns
time = data[:, 0]
concentration = data[:, 1:]

# Function to plot a hexagonal lattice with color-coded concentrations
def plot_hexagonal_lattice(ax, concentration):
    # Hexagonal lattice parameters
    n_rows, n_cols = 20, 20
    radius = 0.5
    orientation = np.sqrt(3) / 2.0
    hex_width = radius * 2
    hex_height = radius * 2 * orientation

    for row in range(n_rows):
        for col in range(n_cols):
            x = col * (3/2) * radius
            y = row * orientation * 2 * radius
            if col % 2 == 1:  # Shift every other row to the right
                y += orientation * radius

            concentration_val = concentration[row * n_cols + col]

            hexagon = RegularPolygon((x, y), numVertices=6, radius=radius, orientation=np.radians(30), edgecolor='black', linewidth=0.5)
            color = plt.cm.viridis(concentration_val)  # Use a colormap to set cell color based on concentration
            hexagon.set_facecolor(color)
            ax.add_patch(hexagon)

    # Calculate the bounds for the x and y axes
    x_margin, y_margin = 0.5, 0.5
    x_min, x_max = -x_margin, (n_cols - 1) * (3/2) * radius + x_margin
    y_min, y_max = -y_margin, (n_rows - 1) * orientation * 2 * radius + orientation * radius + y_margin

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # Set the aspect ratio to equal so that the hexagons are not distorted
    ax.set_aspect('equal', 'box')

# Create a list to store images for the GIF
images = []

# Set the interval for creating frames (every 5 rows)
frame_interval = 100

# Set the time range you want to include in the GIF (from t0 to t1)
t0 = -3
t1 = 200

# Create the colorbar
fig, ax = plt.subplots(figsize=(6, 6))
cax = fig.add_axes([0.87, 0.25, 0.03, 0.5])
cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=0, vmax=2.5)), cax=cax)
cbar.set_label('Concentration', fontsize=12)

# Loop through each time instant within the specified range and create a plot for the hexagonal lattice
for t in range(0, len(time), frame_interval):
    if time[t] > t1:
        break  # Stop when reaching t1

    if t0 <= time[t] <= t1:
        fig, ax = plt.subplots(figsize=(6, 6))
        # Text annotation for the time with no decimal precision
        time_text = ax.text(1.0, -0.05, f'Tiempo: {int(time[t])} u.a.', transform=ax.transAxes, ha='right', fontsize=12)
        plot_hexagonal_lattice(ax, concentration[t])
        plt.axis('off')  # Hide axis for cleaner visualization

        # Save the current plot as a temporary image
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmpfile:
            fig.savefig(tmpfile.name)
            images.append(Image.open(tmpfile.name))

        # Close the current plot for the next iteration
        plt.close()

# Save the images as a GIF
images[0].save("red.gif", save_all=True, append_images=images[1:], duration=200, loop=0)

# Clean up - delete the temporary images
for image in images:
    os.remove(image.filename)
