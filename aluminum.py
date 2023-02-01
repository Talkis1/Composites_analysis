


"""
calculate the thickness and weight needed for aluminum construction

safety factor=1.5

"""


density=2800 #kg/m^3



N_M=[50000/.15,50000,50000/.15,0,0,0] #kN/m

al_y=242e6 #Pa aluminum yield strength, value taken from textbook

Nx=N_M[0]
Ny=N_M[1]
Nxy=N_M[2]

#principal stresses
pr1=(((Nx-Ny)/2)**2+Nxy**2)**.5+((Nx+Ny)/2)
pr2=((Nx+Ny)/2)-(((Nx-Ny)/2)**2+Nxy**2)**.5

sf=1.5

vonmises=(((pr1-pr2)**2+pr2**2+pr1**2)**.5)/(((2**.5)*al_y)/sf)


thick_aluminum=vonmises*10**3

al_weight=vonmises*.15*density

print('Aluminum thickness=',thick_aluminum,'mm','\n'+'Aluminum weight=',al_weight,'kg')


