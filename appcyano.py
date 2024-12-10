#!/usr/bin/env python
# coding: utf-8

# In[3]:

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
a=9.33
b=1.67
c=4.17e-05
d=4.07e-03
e=2.03e-03
f=1.67e-02
P0=0.0005
C0=0.05
h = 0.005*b
dx = h

# Define the differential equations
def dCad(t, C, P, H, I, L):
    return (C * P - C - H * C * C / P)

def dPad(t, C, P, H, I, L):
    return (I - C * P - L * P)

def run_adsimulation(t, C0, P0, args):
    H, I, L = args
    C = np.zeros(len(t))
    P = np.zeros(len(t))
    C[0] = C0
    P[0] = P0

    for i in range(len(t) - 1):
        if P[i]<=1e-8/(b/a):
            C[i]=0
        if C[i]<=0 and P[i]>=1e-8/(b/a):
            C[i]=1e-7/(b/e)
        k11 = dCad(t[i], C[i], P[i], H, I, L)
        k21 = dPad(t[i], C[i], P[i], H, I, L)

        k12 = dCad(t[i] + 0.5 * dx, C[i] + 0.5 * k11 * dx, P[i] + 0.5 * k21 * dx, H, I, L)
        k22 = dPad(t[i] + 0.5 * dx, C[i] + 0.5 * k11 * dx, P[i] + 0.5 * k21 * dx, H, I, L)

        k13 = dCad(t[i] + 0.5 * dx, C[i] + 0.5 * k12 * dx, P[i] + 0.5 * k22 * dx, H, I, L)
        k23 = dPad(t[i] + 0.5 * dx, C[i] + 0.5 * k12 * dx, P[i] + 0.5 * k22 * dx, H, I, L)

        k14 = dCad(t[i] + dx, C[i] + k13 * dx, P[i] + k23 * dx, H, I, L)
        k24 = dPad(t[i] + dx, C[i] + k13 * dx, P[i] + k23 * dx, H, I, L)

        C[i + 1] = C[i] + (1 / 6) * (k11 + 2 * k12 + 2 * k13 + k14) * dx
        P[i + 1] = P[i] + (1 / 6) * (k21 + 2 * k22 + 2 * k23 + k24) * dx

    return C, P

# Streamlit UI
st.set_page_config(layout="wide", page_title="C-P Dynamics Simulation")
st.title("Cyanobacteria-Phosphorus Dynamics")
st.latex(r"\frac{dC}{dt}=CP-C-h\frac{C{^2}}{P}")
st.latex(r"\frac{dP}{dt}=I-CP-lP")
col1, spacer, col2 = st.columns([1, 0.2, 2])
with col1:
    t = st.slider(r"$t$", 0, 1000, 360, 1)
    H = st.slider(r"$h$", 0.00, 10.0, 0.11, 0.001, format="%.5f")
    I = st.slider(r"$I$", 0.00, 1.0, 0.014, 0.001, format="%.5f")
    L = st.slider(r"$l$", 0.00, 1.0, 0.005, 0.001, format="%.5f")
    C0 = st.number_input("Initial Cyanobacteria concentration (g/L)", min_value=0.00001,value=0.05,step=0.00001,format="%.5f")/(b/e)
    
    P0 = st.number_input("Initial Phosphorus concentration (mg/L)",min_value=0.00001,value=0.005,step=0.00001,format="%.5f")/(b/a)
# Time setup
tad = np.arange(0,t*b+h*b, h*b)

# Run simulation
Cad, Pad = run_adsimulation(tad, C0, P0, (H, I, L))
C=Cad*(b/e)
P=Pad*(b/a)
t = tad/b
# Plot results
with col2:
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()
    sns.lineplot(x=t, y=C, ax=ax1, color="g", label=None)
    sns.lineplot(x=t, y=P, ax=ax2, color="r", label=None)
    ax1.set_xlabel("Time (d)", fontsize=12)
    ax1.set_ylabel("Cyanobacteria concentration (g/L)", color="g", fontsize=12)
    ax2.set_ylabel("Phoshorus concentration (mg/L)", color="r", fontsize=12)
    ax1_lines = [plt.Line2D([0], [0], color="g", lw=2, label="C")]
    ax2_lines = [plt.Line2D([0], [0], color="r", lw=2, label="P")]
    ax1.legend(handles=ax1_lines + ax2_lines, loc="upper right", title="Legend")
    ax1.grid(True)
    ax1.set_title("Cyanobacteria-Phosphorus dynamics")
    st.pyplot(fig)
    
    fig2, ax3 = plt.subplots(figsize=(10, 6))
    ax3.plot(P, C)
    ax3.set_xlabel("P")
    ax3.set_ylabel("C")
    ax3.set_title("Phase portrait")
    st.pyplot(fig2)

