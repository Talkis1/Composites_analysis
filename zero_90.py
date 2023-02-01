# [0m/90n]s

import numpy as np
from collections import Counter

"""
Created on Mon Dec  5 16:00:38 2022

alpha=2*m/N
beta=2*n/N
gamma=4*p/N


N=2*(m+n+2*p)

#allowable laminate thickness = 

design safety factor=1.5


x=1 #m

y=.15 #m

#report needs to be ten pages long including figures and references using arial 11, 1.2 line spacing, and 1" margins

#design needs to include alternative designs based on different lamina's stacking sequence and justify design choices: follow first ply
failure

#discuss and compare different designs and suggest most suitable based on weight and price: discuss suggestions. compare results with
design based on 2024-T3 aluminum alloy

#will use S-glass/epoxy, Kevlar (aramid49) epoxy, and carbon (IM7)/epoxy S-glass=$10/kg, Kevlar=$35/kg, carbon=$60/kg

"""

'''
[0m/90n]s
'''



def qbar(E1,v12,E2,G12,theta):
    v21=E2*v12/E1
    q11=E1/(1-v12*v21)
    q12=(v12*E2)/(1-v12*v21)
    q22=E2/(1-v12*v21)
    q66=G12
    Q=np.array([[q11,q12,0],[q12,q22,0],[0,0,q66]])
    theta=np.radians(theta)
    m=np.cos(theta)
    n=np.sin(theta)
    
    T=np.array([[m**2,n**2,2*m*n],[n**2,m**2,-2*m*n],[-m*n,m*n,m**2-n**2]])
    TinvTrans=np.linalg.inv(np.transpose(T))
    Qbar=np.matmul(np.linalg.inv(T),Q)
    Qbar=np.matmul(Qbar,TinvTrans)
    
    return Qbar



def tsai_wu(s1t,s1c,s2t,s2c,t12f,s1,s2,t6):
    F1=(1/s1t)-(1/s1c)
    F2=(1/s2t)-(1/s2c)
    F66=(1/t12f**2)
    F11=1/(s1t*s1c)
    F22=1/(s2t*s2c)
    F12=-.5*((F11*F22)**.5)
    # print(F11,F22)
    a=F11*s1**2+F22*s2**2+F66*t6**2+2*F12*s1*s2
    b=F1*s1+F2*s2
    # print(a,b)
    Sfa=(-b+(b**2+4*a)**.5)/(2*a)
    Sfr=(-b-(b**2+4*a)**.5)/(2*a)
    return Sfa,Sfr

def Aa_mats(angles,k,E1,E2,G12,v12):
    A11=[]
    A12=[]
    A16=[]
    A21=[]
    A22=[]
    A26=[]
    A61=[]
    A62=[]
    A66=[]
    B11=[]
    B12=[]
    B16=[]
    B21=[]
    B22=[]
    B26=[]
    B61=[]
    B62=[]
    B66=[]
    D11=[]
    D12=[]
    D16=[]
    D21=[]
    D22=[]
    D26=[]
    D61=[]
    D62=[]
    D66=[]

    n=0
    for theta in angles:
        
        Qbar=qbar(E1,v12,E2,G12,theta)

        
        
        z=k[n+1]-k[n]

        A11.append(Qbar[0][0]*z)
        A22.append(Qbar[1][1]*z)
        A66.append(Qbar[2][2]*z)
        A12.append(Qbar[0][1]*z)
        A16.append(Qbar[0][2]*z)
        A21.append(Qbar[1][0]*z)
        A26.append(Qbar[1][2]*z)
        A61.append(Qbar[2][0]*z)
        A62.append(Qbar[2][1]*z)

        
        B11.append(.5*Qbar[0][0]*(k[n]**2-k[n+1]**2))
        B22.append(.5*Qbar[1][1]*(k[n]**2-k[n+1]**2))
        B66.append(.5*Qbar[2][2]*(k[n]**2-k[n+1]**2))
        B12.append(.5*Qbar[0][1]*(k[n]**2-k[n+1]**2))
        B16.append(.5*Qbar[0][2]*(k[n]**2-k[n+1]**2))
        B21.append(.5*Qbar[1][0]*(k[n]**2-k[n+1]**2))
        B26.append(.5*Qbar[1][2]*(k[n]**2-k[n+1]**2))
        B61.append(.5*Qbar[2][0]*(k[n]**2-k[n+1]**2))
        B62.append(.5*Qbar[2][1]*(k[n]**2-k[n+1]**2))
        
        D11.append((1/3)*Qbar[0][0]*(k[n+1]**3-k[n]**3))
        D22.append((1/3)*Qbar[1][1]*(k[n+1]**3-k[n]**3))
        D66.append((1/3)*Qbar[2][2]*(k[n+1]**3-k[n]**3))
        D12.append((1/3)*Qbar[0][1]*(k[n+1]**3-k[n]**3))
        D16.append((1/3)*Qbar[0][2]*(k[n+1]**3-k[n]**3))
        D21.append((1/3)*Qbar[1][0]*(k[n+1]**3-k[n]**3))
        D26.append((1/3)*Qbar[1][2]*(k[n+1]**3-k[n]**3))
        D61.append((1/3)*Qbar[2][0]*(k[n+1]**3-k[n]**3))
        D62.append((1/3)*Qbar[2][1]*(k[n+1]**3-k[n]**3))
        n+=1
    ABD=np.array([[sum(A11),sum(A12),sum(A16),sum(B11),sum(B12),sum(B16)],
    [sum(A21),sum(A22),sum(A26),sum(B21),sum(B22),sum(B26)],
    [sum(A61),sum(A62),sum(A66),sum(B61),sum(B62),sum(B66)],
    [sum(B11),sum(B12),sum(B16),sum(D11),sum(D12),sum(D16)],
    [sum(B21),sum(B22),sum(B26),sum(D21),sum(D22),sum(D26)],
    [sum(B61),sum(B62),sum(B66),sum(D61),sum(D62),sum(D66)]])
    for i in np.arange(0,len(ABD+1)):
        for b in np.arange(0,len(ABD+1)):
            if ABD[i][b] < 1e-9:
                ABD[i][b]=0
    A=np.array([[sum(A11),sum(A12),sum(A16)],[sum(A21),sum(A22),sum(A26)],
    [sum(A61),sum(A62),sum(A66)]])
    B=np.array([[sum(B11),sum(B12),sum(B16)],[sum(B21),sum(B22),sum(B26)],
    [sum(B61),sum(B62),sum(B66)]])
    D=np.array([[sum(D11),sum(D12),sum(D16)],[sum(D21),sum(D22),sum(D26)],
    [sum(D61),sum(D62),sum(D66)]])
    Bstar=-np.matmul(np.linalg.inv(A),B)
    Cstar=np.matmul(B,np.linalg.inv(A))
    Dstar=D-np.matmul(Cstar,B)
    a=np.linalg.inv(A)-np.matmul(np.matmul(Bstar,np.linalg.inv(Dstar)),Cstar)
    b=np.matmul(Bstar,np.linalg.inv(Dstar))
    c=-1*np.matmul(np.linalg.inv(Dstar),Cstar)
    d=np.linalg.inv(Dstar)
    # abcd=np.array([[a[0][0],a[0][1],a[0][2],b[0][0],b[0][1],b[0][2]],[a[1][0],a[1][1],a[1][2],b[1][0],b[1]
    # [1],b[1][2]],[a[2][0],a[2][1],a[2][2],b[2][0],b[2][1],b[2][2]],
    #   [c[0][0],c[0][1],c[0][2],d[0][0],d[0][1],d[0][2]],[c[1][0],c[1][1],c[1][2],d[1][0],d[1][1],d[1]
    # [2]],[c[2][0],c[2][1],c[2][2],d[2][0],d[2][1],d[2][2]]])
    abcd=np.linalg.inv(ABD)
    return ABD,abcd


#[S-Glass/epoxy (sg), Kevlar(aramid49)/epoxy (K), carbon(IM7)/epoxy (IM7), 2024-T3 AA (AA)]
E1=[45e9,80e9,190e9,73e9] #pascal
E2=[ 11e9, 5.5e9, 9.9e9, 73e9]
G12=[ 4.5e9, 2.2e9, 7.8e9, 26.6e9]
v12=[ .29, .34, .35, .33]
density=[ 2000, 1380, 1610, 2800] #kg/m**3
F1t=[ 1725e6, 1400e6, 3250e6, 414e6] #pascal
F2t=[ 49e6, 30e6, 62e6, 414e6]
F6=[ 70e6, 49e6, 75e6, 248e6]
F1c=[ 690e6, 335e6, 1590e6, 414e6]
F2c=[ 158e6, 158e6, 200e6, 414e6]

ang_m=0
ang_n=90
# k1=[]
th=[.000165, .000127,.000127]  #m
# for i in
angles=[0, 90, 90, 0]
k=[-th[0]*2,-th[0]*1,0,th[0],th[0]*2]
matrices=Aa_mats(angles, k,E1[0],E2[0],G12[0],v12[0])
ABD1=matrices[0]
abcd1=matrices[1]

N_M=[50000/.15,50000,50000/.15,0,0,0]
# N_M=[460000,920000, 228000,0,0,0]
N_M=np.transpose(N_M)
# it=0
angles_dic={}
for i in np.arange(0,3):
    if i==0:
        typ='sglass'
    if i==1:
        typ='K'
    if i==2:
        typ='IM7'
    
    angles_dic={}
    
    
    for m_ang in np.arange(1,40):
        
        # break
        for n_ang in np.arange(1,40):
                
            # break
                
            angles1=[0 for it in np.arange(0,m_ang)]
            angles2=[90 for it in np.arange(0,n_ang)]
            
    
            stackup=[]
            stackup.extend(angles1)
            stackup.extend(angles2)
            stackup.extend(angles2)
            stackup.extend(angles1)
            angles=stackup
            # print(angles)
            
            numlist=[i*.000165 for i in np.arange(0,int(len(angles)/2)+.001)]
            numlist1=[i*-1*.000165 for i in np.arange(int(len(angles)/2),0,-1)]
            numlist2=[i*.000127 for i in np.arange(0,int(len(angles)/2)+.001)]
            numlist3=[i*-1*.000127 for i in np.arange(int(len(angles)/2),0,-1)]
            numlist1.extend(numlist)
            numlist3.extend(numlist2)
            k1=[]
            k1.append(numlist1)
            k1.append(numlist3)
            k1.append(numlist3)
            k=k1[i]
            # print(k)
            matrices=Aa_mats(angles, k,E1[i],E2[i],G12[i],v12[i])
            ABD=matrices[0]
            abcd=matrices[1]

            
            strains=np.matmul(abcd,N_M)
            
            tstrain=[]
            xstrain=[]
            ystrain=[]
            
            for it in range(0,len(k)-1):
                thick=k[it]/2+k[it+1]/2
        
                xstrain.append(strains[3]*thick+strains[0])
                ystrain.append(strains[4]*thick+strains[1])
                tstrain.append(strains[5]*thick+strains[2])
        
        
            gstress=[]
            lstress=[]
            lstrain=[]
            p=0
            for it in angles:
                strainx=xstrain[p]
        
                strainy=ystrain[p]
                straint=tstrain[p]
                strains=np.array([[strainx],[strainy],[straint]])
        
                q=qbar(E1[i],v12[i],E2[i],G12[i],it)
                gstress.append(np.matmul(q,strains))
                theta=np.radians(it)
                m=np.cos(theta)
                n=np.sin(theta)
                T=np.array([[m**2,n**2,2*m*n],[n**2,m**2,-2*m*n],[-m*n,m*n,m**2-n**2]])
                lstress.append(np.matmul(T,np.matmul(q,strains)))
                lsss=np.matmul(T,np.matmul(q,strains))
                s11=1/E1[i]
                s12=-v12[i]/E1[i]
                s22=1/E2[i]
                s66=1/G12[i]
                S=np.array([[s11,s12,0],[s12,s22,0],[0,0,s66]])
                lstrain.append(np.matmul(S,lsss))
                p+=1
        
            tsai1=[]
            tsai2=[]
            for it in lstress:
                s1=it[0]
                s2=it[1]
                t12=it[2]
        
                tsai=tsai_wu(F1t[i], F1c[i], F2t[i], F2c[i], F6[i], s1, s2, t12)
        
                tsai1.append(tsai[0][0])
                tsai2.append(tsai[1][0])
        
            tsai1=min(tsai1,key=abs)
            tsai2=min(tsai2,key=abs)
            if tsai1<0:
                tsai1=tsai1*-1
            if tsai2<0:
                tsai2=tsai2*-1
            tsai3=[tsai1,tsai2]
            tsai4=min(tsai3,key=abs)
        
            m_initial=1.5/tsai4
            
            if tsai4>=1.5:
                total_thickness=2*(m_ang+n_ang)*th[i]
                weight=.15*total_thickness*2*density[i]
                
                angles_dic[typ]=[m_ang,n_ang,tsai4,m_ang+n_ang,total_thickness,.15*total_thickness*density[i]]
                
                break
        print(angles_dic)










