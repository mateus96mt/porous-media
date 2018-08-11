import numpy as np
import aux
import os

def passo1(c, nx, ny, domx, domy, vel, D11, D12, D21, D22, alpha_mol, alpha_l, alpha_t, phi, h, dt, p, q, t, Tinj):
    
    da = np.zeros(nx+1) #diagonal abaixo da principal
    db = np.zeros(nx+1) #diagonal principal
    dc = np.zeros(nx+1) #diagonal acima da principal
    fd = np.zeros(nx+1) #diagonal principal
    
    for j in range(ny+1):
        
        for i in range(nx+1):
            
            di = D11[i+1,j] if i<nx else D11[i,j] 
            esq = D11[i-1,j] if i<nx else D11[i,j]
            
            
            D11pi12 = (D11[i,j] + di)/2.0
            D11li12 = (D11[i,j] + esq)/2.0
            
            da[i] = -D11li12/(h**2)
            db[i] = (D11pi12 + D11li12)/(h**2) + 2*phi/dt
            dc[i] = -D11pi12/(h**2)
            
            fd[i] = aux.dxdy(i, j, c, h, q, D12, nx, ny) 
            fd[i] += aux.dydx(i, j, c, h, q, D21, nx, ny) 
            fd[i] += aux.dydy(i, j, c, h, q, D22, nx, ny) 
            fd[i] -= aux.upwy(i, j, c, vel, h, q, nx, ny)
            fd[i] += c[i,j,q]*(2*phi/dt)
            
            if vel[i,j,0]>=0:
                da[i] += -vel[i,j,0]/h 
                db[i] += vel[i,j,0]/h
            else:
                dc[i] += vel[i,j,0]/h
                db[i] += -vel[i,j,0]/h
        
            if i==0:
                fd[i] += -da[i]*c[i,j,q]
            if i==nx:
                fd[i] += -dc[i]*c[i,j,q]
        
        '''injeção no ponto 0,0, forcando ser sempre 1 nesse ponto'''
        if j==0:
            db[0] = 1.0
            dc[0] = 0.0
            fd[0] = 1.0
            
            if t>Tinj:
                fd[0] = 0.0
        
        c[:, j, p] = aux.TDMAsolver(da[1:], db, dc[:-1], fd)

def passo2(c, nx, ny, domx, domy, vel, D11, D12, D21, D22, alpha_mol, alpha_l, alpha_t, phi, h, dt, p, q, t, Tinj):
    
    da = np.zeros(ny+1) #diagonal abaixo da principal
    db = np.zeros(ny+1) #diagonal principal
    dc = np.zeros(ny+1) #diagonal acima da principal
    fd = np.zeros(ny+1) #diagonal principal
    
    for i in range(nx+1):
        
        for j in range(ny+1):
            
            cima = D22[i,j+1] if j<nx else D22[i,j] 
            baixo = D22[i,j-1] if j<nx else D22[i,j] 
            
            D22pj12 = (D22[i,j] + cima)/2.0
            D22lj12 = (D22[i,j] + baixo)/2.0
            
            da[j] = -D22lj12/(h**2)
            db[j] = (D22pj12 + D22lj12)/(h**2) + 2*phi/dt
            dc[j] = -D22pj12/(h**2)
            
            fd[j] = aux.dxdy(i, j, c, h, q, D12, nx, ny) 
            fd[j] += aux.dydx(i, j, c, h, q, D21, nx, ny) 
            fd[j] += aux.dxdx(i, j, c, h, q, D11, nx, ny) 
            fd[j] -= aux.upwx(i, j, c, vel, h, q, nx, ny)
            fd[j] += c[i,j,q]*(2*phi/dt)
            
            if vel[i,j,1]>=0:
                da[j] += -vel[i,j,1]/h 
                db[j] += vel[i,j,1]/h
            else:
                dc[j] += vel[i,j,1]/h
                db[j] += -vel[i,j,1]/h
        
            if j==0:
                fd[j] += -da[j]*c[i,j,q]
            if j==ny:
                fd[j] += -dc[j]*c[i,j,q]
        
        '''injeção no ponto 0,0, forcando ser sempre 1 nesse ponto'''
        if i==0 and t<=Tinj:
            db[0] = 1.0
            dc[0] = 0.0
            fd[0] = 1.0
            
            if t>Tinj:
                fd[0] = 0.0
        
        c[i, :, p] = aux.TDMAsolver(da[1:], db, dc[:-1], fd)

      
def ADI_transporte(c, nx, ny, nt, domx, domy, vel, D11, D12, D21, D22, alpha_mol, alpha_l, alpha_t, phi, h, dt, Tinj):
    
    p, q = 1, 0
    
    aux.gera_vtk("dados/t0.vtk", c, nx, ny, p)
    
    for t in range(nt):
        
        p, q = q, p
        passo1(c, nx, ny, domx, domy, vel, D11, D12, D21, D22, alpha_mol, alpha_l, alpha_t, phi, h, dt, p, q, t, Tinj)
        
        p, q = q, p
        passo2(c, nx, ny, domx, domy, vel, D11, D12, D21, D22, alpha_mol, alpha_l, alpha_t, phi, h, dt, p, q, t, Tinj)

        aux.gera_vtk("dados/t" + str(t+1) + ".vtk", c, nx, ny, p)

def main():
    #os.system("gfortran -o ../darcy/main ../darcy/darcy-1m.f")
    #os.system("cd ../darcy;./main -> vel.dat")
    
#    1.8540746773
    domx = [0, 1000] #dominio em x
    domy = [0, 1000] #dominio em y
    domt = [0, 2000] #dominio no tempo
    
    Tinj = domt[1]
    
    nx = 40
    ny = 40
    
#    vel = aux.leCampoVelocidade("../darcy/vel.dat", nx, ny, 0)
    vel = aux.leCampoVelocidade("../darcy/velocity.dat", nx, ny, 1)
    vel[vel<1e-10]=0 #joga pra zero valores muito pequenos
    norm_vel = np.sqrt(vel[:,:,0]**2 + vel[:,:,1]**2)
    norm_vel[norm_vel<1e-10] = 1
    
    hx = (domx[1]-domx[0])/nx
    hy = (domy[1]-domy[0])/ny
    h = hx
    
    dt = 2.5
    
    nt = int((domt[1] - domt[0])/dt)
    
    alpha_mol = 0.0 #difusao molecular
    alpha_l = 1.0 #dispersao longitudinal
    alpha_t = 0.0 #dispersao transversal
    phi = 0.1 #porosidade
    
    #tensor de dispersao
    D11 = np.zeros((nx+1, ny+1))
    D12 = np.zeros((nx+1, ny+1))
    D22 = np.zeros((nx+1, ny+1))
    D11[:,:] = alpha_mol + (alpha_l*vel[:,:,0]**2 + alpha_t*vel[:,:,1]**2)/norm_vel[:,:]
    D22[:,:] = alpha_mol + (alpha_l*vel[:,:,1]**2 + alpha_t*vel[:,:,0]**2)/norm_vel[:,:]
    D12[:,:] = (alpha_l*vel[:,:,0]*vel[:,:,1] - alpha_t*vel[:,:,0]*vel[:,:,1])/norm_vel[:,:]
    D21 = D12
    c = np.zeros((nx+1, ny+1, 2))
    
    
    """teste gerar mesmos resultados DCC190"""
#    D11 = np.full((nx+1, ny+1), 1e-5)
#    D22 = np.full((nx+1, ny+1), 1e-5)
#    D12 = np.zeros((nx+1, ny+1))
#    D21 = D12
#    phi = 1.0
#
#    aux.geraVetores(domx[0], domx[1], domy[0], domy[1], vel)
#    print(vel[0,0,0])
#    print(vel[0,0,1])
#    print("h=", hx, "   dt=", dt, "   nx =", nx, "   ny=", ny, "   nt=", nt)
    """teste gerar mesmos resultados DCC190"""
    
    
    ADI_transporte(c, nx, ny, nt, domx, domy, vel, D11, D12, D21, D22, alpha_mol, alpha_l, alpha_t, phi, h, dt, Tinj)

os.system("rm -r dados")
main()




#nx = 40
#ny = 40
#
##    vel = aux.leCampoVelocidade("../darcy/vel.dat", nx, ny, 0)
#vel = aux.leCampoVelocidade("../darcy/velocity.dat", nx, ny, 1)