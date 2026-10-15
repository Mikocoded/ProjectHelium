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

x = df['x'].to_numpy()
y = df['y'].to_numpy()
z = df['z'].to_numpy()

points = np.column_stack((
    x,
    y,  
    z
))

prob = df['prob'].to_numpy()

if deci1 == 0:
    deci1_1 = int(input("FOr Scatter plot type 0, for Volume plot type 1: "))

    if deci1_1 == 0:
        prob_norm = df['prob'].to_numpy() / df['prob'].max()

        y = points[:, 1]
        y_min, y_max = y.min(), y.max()
        
        pl = pv.Plotter()
        pl.set_background('white')
        def slicer(val) -> None:
            mask = y <= val
            sub = pv.PolyData(points[mask])
            sub['Probability'] = prob_norm[mask]
            pl. add_mesh(
                sub,
                scalars='Probability',
                cmap='inferno',
                clim = [0,1],
                opacity='linear',
                scalar_bar_args={'title': 'Probability', 'color': 'black'},
                name='prob_cloud',
            )
        pl.add_text(f'Helium single-electron probability | Energy = {energia} Hartree',
                    font_size=10, color='white')
        pl.add_slider_widget(slicer, [y_min, y_max], value=y_max, title='Slicer', interaction_event='end')
        log("Rendering")
        pl.show()
    elif deci1_1 == 1:
        lower = np.percentile(points, 0.5, axis=0).astype(np.float32)
        upper = np.percentile(points, 99.5, axis=0).astype(np.float32)

        n = 400

        hist, edges = np.histogramdd(points, bins=n, range=list(zip(lower, upper)), weights=prob.astype(np.float16))
        hist = hist.astype(np.float32, copy=False)
        hist = gaussian_filter(hist, sigma=3.0)

        nonzero = hist[hist > 0]
        vmax = np.percentile(nonzero, 99.5).astype(np.float16)
        del nonzero

        hist /= (vmax * 10)
        np.log1p(hist, out=hist)
        hist /= hist.max()

        grid = pv.ImageData()
        grid.dimensions = np.array(hist.shape)
        grid.origin = lower
        grid.spacing = (upper - lower) / (np.array(hist.shape) - 1)
        grid.point_data['Probability'] = hist.flatten(order='F')
        del hist

        nx, ny, nz = grid.dimensions

        y_profile = grid.point_data['Probability'].reshape((nx, ny, nz), order='F').sum(axis=(0, 2))
        threshold = y_profile.max() * 0.5
        active = np.where(y_profile > threshold)[0]
        y_lo, y_hi = int(active.min()), int(active.max())

        pl = pv.Plotter()
        pl.set_background('white')
        def slicer(value):
            cut = int(value)
            sub = grid.extract_subset([0, nx - 1, 0, cut, 0, nz - 1])
            pl.add_volume(
                sub,
                scalars='Probability',
                cmap='inferno',
                clim=[0,1],
                opacity='linear',
                shade=True,
                scalar_bar_args={'title': 'Probability', 'color': 'black'},
                name='prob_volume',
                )
        pl.add_slider_widget(slicer, [0, ny - 1], value=ny - 1, title='Slicer', interaction_event='end')
        pl.add_text(f'Helium single-electron probability | Energy = {energia} Hartree',
            font_size=10, color='black')
        pl.camera.zoom(10.0)
            # pl.screenshot('helium.png')
        pl.show()
            # log("Saved helium.png")
        



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
