import matplotlib.pyplot as plt
import numpy as np

# X and Y axis edges
y_values = [20, 30, 60, 120]
x_values = [0, 1, 1.8, 2.5]

# Data to be displayed in each quadricule (row-major order)
data = [
    0.66378746, 0.63024427, 0.587,
    0.74141379, 0.7083548, 0.666, 
    0.78504988, 0.74889387, 0.695   #b
   
]

data = [
    0.094, 0.098, 0.092, 
    0.086, 0.097, 0.095, #c
    0.097, 0.103, 0.101
]

data = [
    0.005, 0.006, 0.005, 
    0.004, 0.005, 0.005, 
    0.004, 0.005, 0.006 #light
    
]

# Reshape the data into a 4x4 grid
data_grid = np.array(data).reshape(3, 3)

# Create the plot with a custom figure size
fig, ax = plt.subplots(figsize=(8, 6))  # Widen the plot more to make space for narrower x-values



# Plot the grid with the data, using the pcolormesh function for custom edges
c = ax.pcolormesh(x_values, y_values, data_grid, cmap="coolwarm", shading='flat')

# Annotate each square with the corresponding value
for i in range(data_grid.shape[0]):
    for j in range(data_grid.shape[1]):
        # Find the x and y positions in the center of each square
        x_pos = (x_values[j] + x_values[j+1]) / 2
        y_pos = (y_values[i] + y_values[i+1]) / 2
        ax.text(x_pos, y_pos, f'{data_grid[i, j]:.3f}', ha='center', va='center', fontsize=14, color='black')

# Add colorbar for better visualization
fig.colorbar(c)

# Set axis labels
ax.set_ylabel(r'$p_{\mathrm{T}}$ (GeV)')
ax.set_xlabel(r'$|\eta|$')

# Set equal aspect ratio so that the grid fits properly
#ax.set_aspect('auto')

# Adjust layout to ensure all elements fit
plt.tight_layout()

# Display the plot
#plt.show()



# Display the plot
plt.savefig('plots_beff/btagSF_medium_l.pdf')
#plt.savefig('plots_beff/btagSF_medium_c.png')
#plt.savefig('prueba.pdf')
