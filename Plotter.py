import pandas as pd
import numpy as np
import plotly.graph_objects as go
import pyvista as pv
import time
from scipy.ndimage import gaussian_filter

deci1 = int(input("For local plot type 0, for website type 1: "))

t0 = time.time()
def log(msg):
    print(f"[{time.time()-t0:5.1f}s] {msg}", flush=True)


log("Loading data...")
df = pd.read_csv('HeProbability.csv')

energia = df['energy'].max()

log(f"Plotting {len(df)} points...")

r = df['r'].to_numpy()
theta = df['theta'].to_numpy()
phi = df['phi'].to_numpy()

points = np.column_stack((
    r * np.sin(theta) * np.cos(phi),
    r * np.sin(theta) * np.sin(phi),
    r * np.cos(theta),
))

prob = df['prob'].to_numpy()

if deci1 == 0:
    deci1_1 = int(input("FOr Scatter plot type 0, for Volume plot type 1: "))

    if deci1_1 == 0:
        prob_norm = df['prob'].to_numpy() / df['prob'].max()
        cloud = pv.PolyData(points)
        cloud['Probability'] = prob_norm

        pl = pv.Plotter()
        pl.set_background('black')
        pl.add_mesh(
            cloud,
            scalars='Probability',
            cmap='inferno',
            clim = [0,1],
            opacity='linear',
            scalar_bar_args={'title': 'Probability', 'color': 'white'},
        )
        pl.add_text(f'Helium single-electron probability | Energy = {energia} Hartree',
                    font_size=10, color='white')
        log("Rendering")
        pl.show()
    elif deci1_1 == 1:
        n = 360
        mins, maxs = points.min(axis=0), points.max(axis=0)
        hist, edges = np.histogramdd(points, bins=n, range=list(zip(mins, maxs)), weights=prob)

        hist = gaussian_filter(hist, sigma=1.5)

        hist_log = np.log1p(hist)    
        hist_log /= hist_log.max()

        nonzero = hist[hist > 0]
        vmax = np.percentile(nonzero, 99.5)
        hist_log = np.log1p(hist / vmax * 10) 
        hist_log /= hist_log.max()

        grid = pv.ImageData()
        grid.dimensions = np.array(hist_log.shape)
        grid.origin = mins
        grid.spacing = (maxs - mins) / (np.array(hist_log.shape) - 1)
        grid.point_data['Probability'] = hist_log.flatten(order='F')

        pl = pv.Plotter()
        pl.set_background('black')
        pl.add_volume(
            grid,
            scalars='Probability',
            cmap='inferno',
            clim=[0,1],
            opacity='sigmoid_6',
            shade=True,
            scalar_bar_args={'title': 'Probability', 'color': 'white'},
        )
        pl.add_text(f'Helium single-electron probability | Energy = {energia} Hartree',
            font_size=10, color='white')
        pl.show()
        



elif deci1 == 1:
    fig = go.Figure(data=[go.Scatter3d(
        x = points[:,0],
        y = points[:,1],
        z = points[:,2],
        mode='markers',
        marker=dict(
            size=2,
            color=df['prob'],
            colorscale='Viridis',
            opacity=0.3,
            colorbar=dict(title="Probability")
        )
    )])
    fig.update_layout(
        title=f'Electron claude (Single electron helium atom probability) | Energy={energia} Hartree',
        scene=dict(
            xaxis_title='X (Bohrs radius)',
            yaxis_title='Y (Bohrs radius)',
            zaxis_title='Z (Bohrs radius)'
        ),
        template='plotly_dark',
    )
    fig.show()
