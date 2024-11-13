import MDAnalysis as md
import nglview as ng
from sys import stdout
from openmm.app import *
import ipywidgets as widgets
from IPython.display import display
import warnings

# Suprimir o aviso específico de depreciação do DCDReader
warnings.filterwarnings("ignore", category=DeprecationWarning, message="DCDReader currently makes independent timesteps")

# Carregar o arquivo PDB
pdb_file_path = 'data/polyGV.pdb'
pdb3_file = PDBFile(pdb_file_path)

### 6. Visualização
sys = md.Universe(pdb_file_path, 'data/polyALA_traj.dcd')
view = ng.show_mdanalysis(sys, gui=True)

# Exibir o widget de visualização
display(view)

# Configuração do Tab e Slider do ipywidgets
tab = widgets.Tab(children=[
    widgets.Box(children=[
        widgets.Box(children=[
            widgets.Box(children=[
                widgets.Label(value='step'),
                widgets.IntSlider(value=1, min=-100)
            ])
        ])
    ])
])

display(tab)
