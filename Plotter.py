import pandas as pd
import numpy as np
import plotly.graph_objects as go

print("Wczytywanie danych...")
df = pd.read_csv('/home/miko/ProjectAtom/Heprobability.csv')

prog = df['prob'].max() * 0.01
df = df[df['prob'] > prog]
energia = df['energy'].max()

print(f"Rysowanie {len(df)} punktów...")

# Konwersja ze współrzędnych sferycznych na kartezjańskie (X, Y, Z)
x = df['r'] * np.sin(df['theta']) * np.cos(df['phi'])
y = df['r'] * np.sin(df['theta']) * np.sin(df['phi'])
z = df['r'] * np.cos(df['theta'])

# Tworzenie trójwymiarowego wykresu punktowego (Scatter3d)
fig = go.Figure(data=[go.Scatter3d(
    x=x,
    y=y,
    z=z,
    mode='markers',
    marker=dict(
        size=2,
        color=df['prob'],
        colorscale='Viridis',
        opacity=0.3,
        colorbar=dict(title="Prawdopodobieństwo")
    )
)])

fig.update_layout(
    title=f'Chmura elektronowa (Prawdopodobieństwo 3D) | Energia={energia} Hartree',
    scene=dict(
        xaxis_title='X (promienie Bohra)',
        yaxis_title='Y (promienie Bohra)',
        zaxis_title='Z (promienie Bohra)'
    ),
    template='plotly_dark',
)


fig.show()