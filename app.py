import streamlit as st
import pyvista as pv
from stpyvista import stpyvista
import cadquery as cq
from cadquery import exporters
from gekko import GEKKO
import math
import pandas as pd

#Set up sidebar and parameters t
st.set_page_config(page_title="L-Unit Optimisation Calculator",page_icon=":computer",layout="wide")

df=pd.DataFrame({'Parameters':[23.0,18.0,20.0,0.5,0.15], 'Units':['kN/m3','kN/m3','degrees','','m']},
               index=['Concrete Density',
                      'Soil Density',
                      'Soil Friction Angle',
                      'Base Friction Coeficient',
                      'Chamfer length'])

with st.sidebar:
    st.header('Infomation')
    st.subheader('About')
    st.write('This calculator has been developed for educational purposes only and any results should be independently verified.')
    st.subheader('Parameters')
    df=st.data_editor(df)
    st.subheader('Planned Development')
    st.write('''
    - Add Generative AI Input
    ''')

    st.subheader('Author')
    st.write('''
             - Tom Cartigny
             - Git Hub: [@TomCarts/Beam/Calculator](https://github.com/TomCarts/Beam_Calculator.git)
             ''')

#Streamlit Page Layout
st.title('L-Unit Optimisation')

st.subheader('Inputs',divider='grey')

col1, col2, col3 = st.columns(3)
with col1:
    h=st.number_input('Retained Height (m)',min_value=0.5,max_value=5.0,value=2.0,step=0.1)
with col2:
    P_a =st.number_input('Accidental Impact Load (kN)',min_value=0.0,max_value=500.0,value=150.0,step=0.1)
with col3:
    q=st.number_input('Surchagre (kN/m2)',min_value=0.0,max_value=100.0,value=10.0,step=0.1)

#Extract parameters from sidebar dataframe
p_c=df.iloc[0,0]
p_s=df.iloc[1,0]
phi=df.iloc[2,0]
Fr=df.iloc[3,0]
ch=df.iloc[4,0]
#Soil Horizontal Coeficient 
k_s  = (1-math.sin(phi*3.14/180))

#Gekko opitmisation analysis
m = GEKKO() # Initialize gekko

# Use IPOPT solver (default)
m.options.SOLVER = 3
# Change to parallel linear solver
m.solver_options = ['linear_solver ma97']
# Initialize variables
Bw = m.Var(value=1,lb=1,ub=5) #Base width limits
Bt = m.Var(value=0.3,lb=0.3,ub=0.75) #Base thickness limits
Uh = m.Var(value=0.5,lb=0.5,ub=5) #Upstand height limits
Ut = m.Var(value=0.3,lb=0.3,ub=0.5) #Upstand thickness limits
dc = m.Var(value=0.1,lb=0.05,ub=0.2) #Chamfer diameter limits
L= m.Var(value=2,lb=0.5,ub=5) #Unit Length limits

# Equations
m.Equation(Bt+Uh==h) #L-Unit height
m.Equation((Bt+Uh)/Bw<=1.0) #L-Unit Height / Width Ratio
m.Equation((Bt+Uh)/L<=1.0) #L-unit Height / Length Ratio
m.Equation(((P_a*(Uh+Bt)+(0.5*(Bt+Uh)**2*p_s*L*(k_s)*((Bt+Uh)/3))+0.5*L*(Bt+Uh)**2*q*k_s)/((Bw*Bt*Bw/2+Uh*Ut*Ut/2)*p_c*L+(Uh*(Bw-Ut))*L*p_s*(Ut+(0.5*Uh*(Bw-Ut)))))<=0.95) #Overturning FOS - Stabilsing forces
m.Equation(((P_a+(0.5*(Bt+Uh)**2*p_s*L*k_s)+L*(Bt+Uh)*q*k_s)/(((Bw*Bt+Uh*Ut)*p_c*L+(Uh*(Bw-Ut))*L*p_c)*Fr))<=0.95)                           #Sliding FOS - Stabilising Forces
#m.Equation((Bw*Bt+Uh*Ut+dc**2/2)*25*L<=50)   #Weight Limits
m.Minimize((Bw*Bt+Uh*Ut+dc**2/2)*25*L) # Equation to minimise - L-unit Weight
m.options.IMODE = 3 # Steady state optimization
m.solve(disp=False) # Solve

#Streamlit results layout
st.subheader('Opitimised L-Unit parameters: Minimising Volume / Weight',divider='grey')
st.write('Minimum Weight: ' + str(math.floor(m.options.objfcnval))+'kN')

l_df=pd.DataFrame({'Parameters':[Bw.value[0],Bt.value[0],Uh.value[0],Ut.value[0],dc.value[0],L.value[0]], 'Units':['m','m','m','m','m','m']},
               index=['Base Width',
                      'Base Thickness',
                      'Upstand Height',
                      'Upstand Thickness',
                      'Chamfer width',
                      'Length'])

#Extract and Summarise Calculation values
Bw = Bw.value[0]
Bt = Bt.value[0]
Uh = Uh.value[0]
Ut = Ut.value[0]
dc = dc.value[0]
L= L.value[0]

Overturning_FOS = (((P_a*(Uh+Bt)+(0.5*(Bt+Uh)**2*p_s*L*(k_s)*((Bt+Uh)/3))+0.5*L*(Bt+Uh)**2*q*k_s)/((Bw*Bt*Bw/2+Uh*Ut*Ut/2)*p_c*L+(Uh*(Bw-Ut))*L*p_s*(Ut+(0.5*Uh*(Bw-Ut)))))) #Overturning FOS - Stabilsing forces

Sliding_FOS = (((P_a+(0.5*(Bt+Uh)**2*p_s*L*k_s)+L*(Bt+Uh)*q*k_s)/(((Bw*Bt+Uh*Ut)*p_c*L+(Uh*(Bw-Ut))*L*p_c)*Fr)))                           #Sliding FOS - Stabilising Forces

# Define CAD Query Geometry f
pts = [
    (0, 0),
    (Bw,0),
    (Bw,Bt),
    (0+Ut+ch,Bt),
    (0+Ut,Bt+ch),   
    (0+Ut,Bt+Uh),
    (0,Bt+Uh),
    (0,0),
]

#Build Model
result = cq.Workplane("XZ").polyline(pts).close().extrude(L)

# Export to DXF
exporters.export(result, "mesh_v1.stl")
exporters.exportDXF(result, "L-unit.dxf")

## Initialize a plotter object
plotter = pv.Plotter(window_size=[800,400])

#Read CAD quer output
mesh=pv.read('mesh_v1.stl')

## Add mesh to the plotter
plotter.add_mesh(mesh, 
                    #scalars='my_scalar', cmap='bwr'
                )

## Final touches
plotter.view_isometric()
#plotter.add_scalar_bar()
plotter.background_color = 'black'

## Plot Mesh
col1, col2 = st.columns([3,1])

with col1:
    stpyvista(plotter
              #panel_kwargs=dict(
                  #orientation_widget=True,
                  #interactive_orientation_widget=True
                  #)
              #key="pv_cube" #Pass a key to avoid re-rendering at each page change
    )
 
with open("L-unit.dxf", "rb") as file:
    dxf_bytes = file.read()

    st.download_button(
        label="Download DXF",
        data=dxf_bytes,
        file_name="L-unit.dxf",
        mime="application/dxf"
    )

with open("mesh_v1.stl", "rb") as file:
    dxf_bytes = file.read()

    st.download_button(
        label="Download stl",
        data=dxf_bytes,
        file_name="mesh_v1.stl",
        mime="application/octet-stream"
    )
    
with col2:
    st.dataframe(l_df)
    
st.subheader('Calculation Summary',divider='grey')
st.write('Overturning Utilisation = ' + str(math.floor(Overturning_FOS*100)) + "% (Unfactored)")
st.write('Sliding Utilisation = ' + str(math.floor(Sliding_FOS*100)) + "% (Unfactored)")
