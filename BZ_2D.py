import io
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

st.title("2D Brillouin Zone Calculator")
st.write("Adjust the crystal structure with the sliders in the sidebar.")

plt.rcParams.update({
    "text.usetex": False,
    "mathtext.fontset": "cm",
    "font.family": "serif",
})

Nx = 100
N = 3

angle_deg = st.sidebar.slider("Lattice Angle (Degrees)", min_value=30, max_value=150, value=90, step=5)

mag_a1 = 1
mag_ratio = st.sidebar.slider("a2/a1", min_value=0.2, max_value=2.0, value=1.0, step=0.1)
mag_a2 = mag_ratio*mag_a1

a1 = np.array([mag_a1,0,0])
a2 = np.array([mag_a2*np.cos(np.radians(angle_deg)),mag_a2*np.sin(np.radians(angle_deg)),0])
z = np.array([0,0,1])

V = np.dot(a1,np.cross(a2,z))
b1 = 2*np.pi*np.cross(a2,z)/V
b2 = 2*np.pi*np.cross(z,a1)/V
kx = np.linspace(-3.1*np.pi/mag_a1,3.1*np.pi/mag_a1,Nx)

fig, ax = plt.subplots(figsize=(4.2, 4))

G_points = []
for h in range(-N,N+1):
    for k in range(-N,N+1):
        G = h*b1 + k*b2
        G_points.append(G)
        normG = np.linalg.norm(G)
        if G[1] != 0:
            ky = (normG**2-2*G[0]*kx) / (2*G[1])
            plt.plot(kx*mag_a1/np.pi,ky*mag_a2/np.pi,'k')
        else:
            xline = np.array([normG**2 / (2*G[0]), normG**2 / (2*G[0])])
            yline = np.array([-3*np.pi/mag_a1, 3*np.pi/mag_a2])
            plt.plot(xline*mag_a1/np.pi,yline*mag_a2/np.pi,'k')

G_pointsArray = np.array(G_points)
plt.plot(G_pointsArray[:,0]*mag_a1/np.pi, G_pointsArray[:,1]*mag_a2/np.pi, 'r.',markersize=10)

ax.set_xlim([-3, 3])
ax.set_ylim([-3, 3])
ax.set_aspect('equal')
ax.set_xlabel(r'$k_x a_1 / \pi$')
ax.set_ylabel(r'$k_y a_2 / \pi$')
ax.set_title(f'Brillouin Zone')
st.pyplot(fig,use_container_width=True)

pdf_buffer = io.BytesIO()
fig.savefig(pdf_buffer, format="pdf",bbox_inches=None)
pdf_buffer.seek(0)
st.download_button(
    label="💾 Save as PDF",
    data=pdf_buffer,
    file_name=f"Brillouin_Zone.pdf",
    mime="application/pdf"
)