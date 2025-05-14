import numpy as np
from matplotlib import pyplot as plt
import sys

def poisson(c_list):
    Phi=np.zeros(N,dtype=float)
    for i in range(1,12):
        ql=0
        qr=0
        for j in range(0,i):
            ql= ql+c_list[j]
        for j in range(i,213):
            qr=qr+c_list[j]
        Phi[i]=Phi[i-1]+(-ql+qr)/epsilon_0*1E-10/epsilon_r1
    for i in range(12,213):
        ql=0
        qr=0
        for j in range(0,i):
            ql= ql+c_list[j]
        for j in range(i,213):
            qr=qr+c_list[j]
        Phi[i]=Phi[i-1]+(-ql+qr)/epsilon_0*1E-10/epsilon_r2
    return Phi

with open("dos", "r") as data:
    lines = data.readlines()
dos = []
for i in range(0,len(lines)):
    dos.append([float(x) for x in np.array(lines[i].split())])
dos=np.array(dos)

def n_e(dos, mu, T):
    kT=T*8.625E-5
    n=0
    de=dos[1][0]-dos[0][0]
    for i in range(0, len(dos)):
        if (dos[i][0]-mu)/kT>60:
            n=n-dos[i][1]*de*0
        elif (dos[i][0]-mu)/kT<-60:
            n=n-dos[i][1]*de
        else:
            n=n-dos[i][1]/(1+np.exp((dos[i][0]-mu)/kT))*de
    return n/4

def VBM(LAO,STO,L):
    vbm=np.zeros(N,dtype=float)
    for i in range(0,12):
        vbm[i]=LAO
    for i in range(12,12+L):
        vbm[i]=STO-(i-12-L)**2/L**2*0.15
    for i in range(12+L,213):
        vbm[i]=STO
    return vbm

def poisson_e(n,delta_ne,mu,dos,T,CBM,VBM,Phi_i,Phi_field):
    num_e=np.zeros(N,dtype=float)
    if ((CBM-Phi_i-Phi_field)[12]-mu)<0:
        #num_e[12]=((CBM-Phi_i-Phi_field)[12]-mu+0.2)/200/1E-10*epsilon_r2*epsilon_0/2
        num_e[12]=((CBM-Phi_i-Phi_field)[12]-mu)/200/1E-10*epsilon_r2*epsilon_0/2
        num_e[N-1]=-num_e[12]
        step=0
        while step<n:
            phi=Phi_i+Phi_field+poisson(num_e)-poisson(num_e)[N-1]
            N_e=np.zeros(N,dtype=float)
            step =step +1
            for i in range(12,60):
                N_e[i]=n_e(dos,mu+phi[i]-CBM[i],T)
            N_e[N-1]=-np.sum(N_e)
            test=np.abs(np.sum(N_e[12]))-np.abs(np.sum(num_e[12]))
            num_e=(N_e-num_e)/10+num_e

            if step==n-1:
                print("unconv")
                break
            if np.abs(test)<delta_ne:
                #print("conv")
                break
    else:
        step=0
        while step<n:
            phi=Phi_i+Phi_field+poisson(num_e)-poisson(num_e)[N-1]
            N_e=np.zeros(N,dtype=float)
            step =step +1
            for i in range(12,60):
                N_e[i]=n_e(dos,mu+phi[i]-CBM[i],T)
            N_e[N-1]=-np.sum(N_e)
            test=np.abs(np.sum(N_e[12]))-np.abs(np.sum(num_e[12]))
            num_e=N_e

            if step==n-1:
                print("unconv")
                break
            if np.abs(test)<delta_ne:
                #print("conv")
                break
    return num_e

def n_i(L1,L2,NL):
    n_i=np.zeros(N,dtype=float)
    n_i[0]=-0.5
    n_i[11]=0.5
    for i in range(0,L1):
        n_i[12-L1+i]=n_i[12-L1+i]-NL/L1
    for i in range(0,L2):
        n_i[12+i]=n_i[12+i]+(NL-N_smear)/L2
    return n_i

def Phi_Field(Volt,L_vac):
    y=Volt/(12*epsilon_r2/200/epsilon_r1+1+L_vac*epsilon_r2/200)
    x=Volt-L_vac*epsilon_r2/200*y
    Phi=np.zeros(N,dtype=float)
    for i in range(0,13):
        Phi[i]=-(x-y)/12*i+x
    for i in range(13,N):
        Phi[i]=-(y)/200*(i-12)+y
    return Phi

def main_new(T,dos,N_i,n_d_s,dx,dt,Phi_field,N_d,mu):
    kT=T*8.625E-5
    if Phi_field[0]<0:
        n_d_step=n_d_s*tes
    else:
        n_d_step=n_d_s
    for i in range(0,n_d_step):
        Phi_d=poisson(N_d+N_i)
        N_e=poisson_e(200,0.0001,mu,dos,T,cbm,vbm,Phi_d,Phi_field)
        Phi=poisson(N_e+N_i+N_d)+Phi_field
        N_d_new=np.zeros(N,dtype=float)
        J_d=np.zeros(N,dtype=float)
        J_e=np.zeros(N,dtype=float)
        for j in range(0,N-1):
            J_d[j]=-kT*(N_d[j+1]-N_d[j])/dx
            if Phi[j+1]-Phi[j]>0:
                J_e[j]=-N_d[j+1]*(Phi[j+1]-Phi[j])/dx
            else:
                J_e[j]=-N_d[j]*(Phi[j+1]-Phi[j])/dx
        J=J_d+J_e
        delta_t=np.zeros(N+1,dtype=float)
        for k in range(1,N-1):
            if (J[k-1]-J[k])==0 or N_d[k]==0:
                delta_t[k]=dt
            elif ((J[k-1]-J[k])/dx)>0:
                delta_t[k]=dt
            else:
                delta_t[k]=-N_d[k]/((J[k-1]-J[k])/dx)
        if (-J[0])==0 or N_d[0]==0:
            delta_t[0]=dt
        elif ((-J[0])/dx)>0:
            delta_t[0]=dt
        else:
            delta_t[0]=-N_d[0]/((-J[0])/dx)
            
        if (J[N-2])==0 or N_d[N-1]==0:
            delta_t[N-1]=dt
        elif ((J[N-2])/dx)>0:
            delta_t[N-1]=dt
        else:
            delta_t[N-1]=-N_d[N-1]/((J[N-2])/dx)
        
        delta_t[N]=dt   
        dt_new=np.min(delta_t)
        #print(dt_new)
        for j in range(1,N-1):
            N_d_new[j]=N_d[j]+(J[j-1]-J[j])/dx*dt_new
        N_d_new[0]=N_d[0]+(-J[0])/dx*dt_new
        N_d_new[N-1]=N_d[N-1]+(J[N-2])/dx*dt_new
        N_d=N_d_new
        #plt.xlim(0,24) 
        #plt.plot(z,N_d)
        #plt.show()
    return N_d, N_e

def write(X):
    with open("test","w") as f:
        for i in range(0,len(X)):
            f.write(str(X[i])+" "+"\n")
    return

epsilon_0 = 8.854E-12
epsilon_r1= 27
epsilon_r2= 10000
N=213
N_smear=0.36
mu=1.73
tes=75
T=3

N_i=n_i(4,28,-0.05)
N_vO_ini=np.zeros(N,dtype=float)
N_vO_ini[0]=N_smear

cbm=VBM(5.8,3.2,24)
vbm=VBM(0.4,0.0,24)
Time = []
N_list=[]


#1. Set the bias voltage and position of the tip
#phi_f=Phi_Field(0,2)
#phi_f=Phi_Field(10,2)
phi_f=Phi_Field(-10,1)

#2. Read the initial concentration distribution of VO
with open("N_d_51000", "r") as data:
    lines = data.readlines()
n_d = []
for i in range(0,len(lines)):
    n_d.append(float(lines[i].split()[0]))
n_d=np.array(n_d)

#3. Run the main iteration equation and set the number of steps
for t in range(0,1000):
    Time.append(t)
    N_d,N_e=main_new(T,dos,N_i,1,1E-10,1E-19,phi_f,n_d,mu)
    a=0
    b=0
    c=0
    d=0
    n_d=N_d
    print((t,np.sum(n_d)))
    for i in range(0,12):
        a=a+N_e[i]
        c=c+N_d[i]
    for i in range(12,N-1):
        b=b+N_e[i]
        d=d+N_d[i]
    N_list.append(" "+str(-a)+" "+str(-b)+" "+str(c/2)+" "+str(d/2)+" "+"\n")

#4. Print the output file
with open("N_d","w") as f:
    for i in range(0,len(n_d)):
        f.write(str(n_d[i]/32*1E16)+" "+"\n")

with open("N_list","w") as f:
    for i in range(0,len(N_list)):
        f.write(str(N_list[i]))

with open("N_e","w") as f:
    for i in range(0,len(N_e)):
        f.write(str(-N_e[i]/16*1E16)+" "+"\n")

z= np.arange(0,N,1)
plt.plot(z,cbm-poisson(N_e+n_d+N_i)+poisson(N_e)[N-1]-phi_f)
plt.plot(z,vbm-poisson(N_e+n_d+N_i)+poisson(N_e)[N-1]-phi_f)
write(vbm-poisson(N_e+n_d+N_i)+poisson(N_e)[N-1]-phi_f+1.2128833102999166)
plt.show()