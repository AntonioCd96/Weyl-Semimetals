#!/usr/bin/env python
# coding: utf-8

# # <center> Hamiltoniano para semimetales de Weyl tipo I --Minimal Models </center>

# El **objetivo general** de este notebook es explorar el Hamiltoniano presentado en el artículo [1] 
# 
# Los conceptos a introducir serán:
# * Hamiltoniano para un semimetal del Weyl 
# * Relación de dispersión generada por este tipo de materiales
# 
# A diferencia del hamiltoniano explorado anteriormente, este econtiene un parametro $\gamma$, que permite estudiar la transición de fase de un semimetal de Weyl tipo I  a uno tipo II. Adicionalmente, este notebook se enfocara sólo en la fase correspondiente al semimemtal de Weyl tipo I ($\gamma$ = 0).
# 
# ---
# <sup>Fuente:  T. M. McCormick, I. Kimchi, and N. Trivedi. Minimal Models for Topological Weyl Semimetals.Phys. Rev. B, 95(7):075133, Feb 2017</sup>

# ## Multiprocesing 
# 

# In[1]:


import multiprocessing as mp
import plotly.graph_objects as go


# In[2]:


get_ipython().run_line_magic('pylab', 'inline')
get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


def EigenV(k):
    k_x,k_y,k_z=k
    E=eigvalsh(HWeyl(k_x,k_y,k_z))
    return E


# In[4]:


res=pi/101 #resolucion

k_xb,k_yb,k_zb=arange(-pi,pi,res),arange(-pi,pi,res),arange(-pi,pi,res)

a  = 1
k_0= pi/2
   #Weyl positions
tx = 1.0/2      
t  = 1.0/2
m  = 2*t
γ  = 0
def HWeyl(k_x,k_y,k_z):   
    HW = array([[γ*(cos(k_x)-cos(k_0))-2*t*sin(k_z),   -(m*(2-cos(k_y)-cos(k_z))+2*tx*(cos(k_x)-cos(k_0)))+2J*t*sin(k_y)],
                [ -(m*(2-cos(k_y)-cos(k_z))+2*tx*(cos(k_x)-cos(k_0)))-2J*t*sin(k_y), γ*(cos(k_x)-cos(k_0))+2*t*sin(k_z)]])
    return HW


# In[5]:


a_d= len(k_xb) #dimension del arreglo
KX,KZ = meshgrid(k_xb,k_zb)
KX    = KX.reshape((a_d*a_d,))
KZ    = KZ.reshape((a_d*a_d,))

k     = column_stack((KX,zeros_like(KX),KZ))


# In[6]:


get_ipython().run_cell_magic('time', '', '\nEk = map(EigenV,k) #función  y los valores que toma\nEk = array(list(Ek))\nprint(Ek)')


# In[7]:


Enm = Ek.T[0].reshape((a_d,a_d)).T#primer T para +/- segundo para X->Z
Enp = Ek.T[1].reshape((a_d,a_d)).T


KX,KZ = meshgrid(k_xb,k_zb)


# In[8]:


DATA = [ go.Surface( z=Enm, x=(KX),y=(KZ),opacity=0.9,  colorbar_x=0.75,colorscale='deep'),
        go.Surface( z=Enp,x=KX,y=KZ,opacity=0.6, colorbar_x=0.9)]


# In[9]:


fig = go.Figure( data=DATA )

fig.update_layout( autosize=False,
                   width = 800, height = 500,
                   margin= dict(l=65, r=50, b=65, t=90),
                   scene = dict(xaxis_title="kx", 
                                yaxis_title="ky", 
                                zaxis_title="E [t]", 
                                xaxis = dict(showbackground=False), 
                                yaxis = dict(showbackground=False),
                                zaxis = dict(showbackground=False)))

fig.show()


# En esta figura se presenta una realcion de dispersion propia de un semimetal de Weyl tipo I. Dicha relación se caracteriza porque la banda de conducción y valencia son simetricas al plano $K_x, K_y$

# In[ ]:




