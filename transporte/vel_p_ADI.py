import numpy as np

#''' compilar o bacalhau, executar e salvar o resultado no arquivo "vel"  '''
#os.system("cd codigo\ modificado;gfortran -o main darcy-1m.f;./main -> vel")

nelh = 64
nelv = 64

vel = open("../darcy/vel")
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

Mp = np.zeros((nelv+1, nelh+1, 2))

for i in range(nelv):
    for j in range(nelh):
        
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
for i in range(nelv+1):
    for j in range(nelh+1):
        texto.append(str(Mp[i, j, 0]) + '  ' + str(Mp[i, j, 1]) + "\n")

vel_field.writelines(texto)
vel_field.close()