import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

def geraVetores(x0, xn, y0, yn, vel_field, cor='g'):
    
    hx, hy = (xn-x0)/(vel_field.shape[0]-1), (yn-y0)/(vel_field.shape[1]-1)
    
    x, y = [x0+hx*i for i in range(vel_field.shape[0])], [y0+hy*i for i in range(vel_field.shape[0])]
        
    X, Y = np.meshgrid(x, y)
    vx = vel_field[:,:,0]#/np.sqrt(vel_field[:,:,0]**2+vel_field[:,:,1]**2)
    vy = vel_field[:,:,1]#/np.sqrt(vel_field[:,:,0]**2+vel_field[:,:,1]**2)
    
    fig = plt.figure()
    plt.quiver(Y[::1, ::1], X[::1, ::1], vx[::1, ::1], vy[::1, ::1], color=cor)
    fig.savefig("vel_field.png", dpi=200)
    plt.show()

'''formato=0 campo de velocidade gerado por elementos finitos
formato=1 campo de velocidade no formato da disciplina DCC190'''
def leCampoVelocidade(nomeArq, nelh, nelv, formato):
    
    '''campo gerado por elementos finitos'''
    if formato==0:
        nomeArq = open(str(nomeArq))
        nomeArq = nomeArq.read().split('\n')
        Me = np.array(nomeArq[:-2]).astype(float)
        
        #agrupa por ponto
        Me = np.array(np.split(Me, Me.shape[0]/2))
        
        #agrupa por elemento
        Me = np.array(np.split(Me, nelh*nelv))
        
        #agrupa por linha
        Me = np.array(np.split(Me, nelv))
        
        #jogar pra zero valores muito pequenos
        Me[abs(Me)<1e-13] = 0.0
        
        '''print("\n\nMe:\n", Me)'''
        
        vel = np.zeros((nelh+1, nelv+1, 2))
        
        for i in range(nelv):
            for j in range(nelh):
                
                vel[i, j, 0] += Me[i, j, 0, 0]
                vel[i, j, 1] += Me[i, j, 0, 1]
                
                vel[i, j+1, 0] += Me[i, j, 1, 0]
                vel[i, j+1, 1] += Me[i, j, 1, 1]
                
                vel[i+1, j+1, 0] += Me[i, j, 2, 0]
                vel[i+1, j+1, 1] += Me[i, j, 2, 1]
                
                vel[i+1, j, 0] += Me[i, j, 3, 0]
                vel[i+1, j, 1] += Me[i, j, 3, 1]
        
        vel[1:-1, 1:-1, :] = vel[1:-1, 1:-1, :]/4.0
        
        vel[1:-1, 0, :] = vel[1:-1, 0, :]/2.0
        vel[1:-1, -1, :] = vel[1:-1, -1, :]/2.0
        
        vel[0, 1:-1, :] = vel[0, 1:-1, :]/2.0
        vel[-1, 1:-1, :] = vel[-1, 1:-1, :]/2.0
    
        return vel

    '''campo gerado no formato da disciplina DCC190'''
    if formato==1:
        result = pd.read_table(str(nomeArq), header = 0, sep = '  ')    
        result = np.array(result)    
        
        
        vel_field = result
        vx = [[ vel_field[i, 0] for i in range(j*(nelh+1), (j+1)*(nelh+1))] for j in range(nelv+1)]
        vy = [[ vel_field[i, 1] for i in range(j*(nelh+1), (j+1)*(nelh+1))] for j in range(nelv+1)]
        vx = np.matrix(vx).transpose()
        vy = np.matrix(vy).transpose()
        
        vel = np.zeros((nelh+1, nelv+1, 2))
        vel[:,:,0] = vx[:,:]
        vel[:,:,1] = vy[:,:]
        
        return vel

def vel_formatoADI(nomeArq, nelh, nelv):

	vel = open(nomeArq)
	vel = vel.read().split('\n')
	Me = np.array(vel[:-2]).astype(float)

	#agrupa por ponto
	Me = np.array(np.split(Me, Me.shape[0]/2))

	#agrupa por elemento
	Me = np.array(np.split(Me, nelh*nelv))

	#agrupa por linha
	Me = np.array(np.split(Me, nelv))

	#jogar pra zero valores muito pequenos
	Me[abs(Me)<1e-13] = 0.0

	'''print("\n\nMe:\n", Me)'''

	Mp = np.zeros((nelh+1, nelv+1, 2))

	for i in range(nelh):
		for j in range(nelv):
			
			Mp[i, j, 0] += Me[i, j, 0, 0]
			Mp[i, j, 1] += Me[i, j, 0, 1]
			
			Mp[i, j+1, 0] += Me[i, j, 1, 0]
			Mp[i, j+1, 1] += Me[i, j, 1, 1]
			
			Mp[i+1, j+1, 0] += Me[i, j, 2, 0]
			Mp[i+1, j+1, 1] += Me[i, j, 2, 1]
			
			Mp[i+1, j, 0] += Me[i, j, 3, 0]
			Mp[i+1, j, 1] += Me[i, j, 3, 1]



	Mp[1:-1, 1:-1, :] = Mp[1:-1, 1:-1, :]/4.0

	Mp[1:-1, 0, :] = Mp[1:-1, 0, :]/2.0
	Mp[1:-1, -1, :] = Mp[1:-1, -1, :]/2.0

	Mp[0, 1:-1, :] = Mp[0, 1:-1, :]/2.0
	Mp[-1, 1:-1, :] = Mp[-1, 1:-1, :]/2.0

	#gera campo de velocidade no formato do codigo ADI
	vel_field = open("vel_field.dat", "w")
	texto = []
	texto.append("vx  vy\n")
	for i in range(nelh+1):
		for j in range(nelv+1):
			texto.append(str(Mp[i, j, 0]) + '  ' + str(Mp[i, j, 1]) + "\n")

	vel_field.writelines(texto)
	vel_field.close()

def TDMAsolver(a, b, c, d):
    nf = len(d) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
        	    
    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc

def gera_vtk(nome, dado, nx, ny, idx):
    print(nome)
    if not os.path.exists(nome.split('/')[0]):
        os.makedirs(nome.split('/')[0])
    
    arqvtk = open(nome, "w")
    texto = []
    
    texto.append("# vtk DataFile Version 3.0\n")
    texto.append("vtk output\n")
    texto.append("ASCII\n")
    texto.append("DATASET RECTILINEAR_GRID\n")
    texto.append("DIMENSIONS " + str(nx+1) + " " + str(ny+1) + " 1\n")
    
    texto.append("X_COORDINATES " + str(nx+1) + " double\n")
    for ix in range(nx+1):
        texto.append(str(ix) + " ")
    texto.append("\n")
    
    texto.append("Y_COORDINATES " + str(ny+1) + " double\n")
    for jy in range(ny+1):
        texto.append(str(jy) + " ")
    texto.append("\n")
    
    texto.append("Z_COORDINATES 1 double\n")
    texto.append("0\n")

    texto.append("POINT_DATA " + str((nx+1)*(ny+1)) + "\n")
    texto.append("FIELD FieldData 1\n")
    texto.append("Concentracao 1 " + str((nx+1)*(ny+1)) + " double\n")
    for jy in range(ny+1):
        for ix in range(nx+1):
                texto.append(str(dado[ix, jy, idx]) + " ")
    texto.append("\n")
    
    arqvtk.writelines(texto)
    arqvtk.close()

def upwx(i, j, c, vel, h, idx, nx, ny):
    if vel[i,j,0]>=0:
        if i==0:
            return 0
        esq = c[i-1,j,idx]
        return vel[i,j,0]*(c[i,j,idx] - esq)/h
    else:
        if i==nx:
            return 0
        di = c[i+1,j,idx]
        return vel[i,j,0]*(di - c[i,j,idx])/h

def upwy(i, j, c, vel, h, idx, nx, ny):
    if vel[i,j,1]>=0:
        if j==0:
            return 0
        esq = c[i,j-1,idx]
        return vel[i,j,1]*(c[i,j,idx] - esq)/h
    else:
        if j==ny:
            return 0
        di = c[i,j+1,idx]
        return vel[i,j,1]*(di - c[i,j,idx])/h

def dxdx(i, j, c, h, idx, D, nx, ny):
    esq = c[i,j,idx]
    di = c[i,j,idx]
    
    Dli12 = D[i,j]
    Dpi12 = D[i,j]
    
    if i<nx:
        di = c[i+1,j,idx]
        Dpi12 = (D[i,j] + D[i+1,j])/2.0
    if i>0:
        esq = c[i-1,j,idx]     
        Dli12 = (D[i,j] + D[i-1,j])/2.0
    
    return (1.0/h**2)*(Dpi12*di - (Dpi12 + Dli12)*c[i,j,idx] + Dli12*esq)

def dxdy(i, j, c, h, idx, D, nx, ny):
    di_cima = c[i,j,idx]
    di_baixo = c[i,j,idx]
    esq_cima = c[i,j,idx]
    esq_baixo =c[i,j,idx]
    Ddi = D[i, j]
    Desq = D[i, j]
    
    if i<nx and j<ny:
        di_cima = c[i+1,j+1,idx]
    if i<nx and j >0:
        di_baixo = c[i+1,j-1,idx]
    if i>0 and j<ny:
        esq_cima = c[i-1,j+1,idx]
    if i>0 and j>0:
        esq_baixo =c[i-1,j-1,idx]

    if i<nx:
        Ddi = D[i+1,j]
    if i>0:
        Desq = D[i-1,j]
    
    return (1.0/(4*(h**2)))*(Ddi*di_cima - Ddi*di_baixo -Desq*esq_cima + Desq*esq_baixo)

def dydx(i, j, c, h, idx, D, nx, ny):
    di_cima = c[i,j,idx]
    di_baixo = c[i,j,idx]
    esq_cima = c[i,j,idx]
    esq_baixo =c[i,j,idx]
    Dcima = D[i, j]
    Dbaixo = D[i, j]
    
    if i<nx and j<ny:
        di_cima = c[i+1,j+1,idx]
    if i<nx and j >0:
        di_baixo = c[i+1,j-1,idx]
    if i>0 and j<ny:
        esq_cima = c[i-1,j+1,idx]
    if i>0 and j>0:
        esq_baixo =c[i-1,j-1,idx]
    
    if j<ny:
        Dcima = D[i,j+1]
    if j>0:
        Dbaixo = D[i,j-1]
    
    return (1.0/(4*(h**2)))*(Dcima*di_cima - Dcima*esq_cima - Dbaixo*di_baixo + Dbaixo*esq_baixo)

def dydy(i, j, c, h, idx, D, nx, ny):
    cima = c[i,j,idx]
    baixo = c[i,j,idx]
    Dpj12 = D[i,j]
    Dlj12 = D[i,j]
    
    if j<ny:
        cima = c[i,j+1,idx]
        Dpj12 = (D[i,j] + D[i,j+1])/2.0
    
    if j>0:
        baixo = c[i,j-1,idx]
        Dlj12 = (D[i,j] + D[i,j-1])/2.0

    return (1.0/h**2)*(Dpj12*cima - (Dpj12 + Dlj12)*c[i,j,idx] + Dlj12*baixo)