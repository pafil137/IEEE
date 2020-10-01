#BIBLI
import csv
import math
import numpy as np

#PLOTS
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

'''
BConfigMatrix = len(_a)*len(_dist)*len(_h)

BConfigMatrix[0] = len(_R)

BConfigMatrix[0][0]
for i in range(len(Blist)):
    for j in range(len(Blist[i])):        
            print("Eihop", Blist[i][j][0])
            print("h", Blist[i][j][1])
            print("d", Blist[i][j][2])
            print("R", Blist[i][j][3])
            print("P0", Blist[i][j][4])
            print("Nb", Blist[i][j][5])
            print("a", Blist[i][j][6])
'''

###########################################################################
##montar Conjunto de treinamento
###########################################################################

#Gerar treino entrada e saida DNN_training NORMAL
def gerarConjuntoTreino(BConfigMatrix, name):

    fileTreino = csv.writer(open("./{0}.csv".format(name), "w"))
    
    fileTreino.writerow(["BEnh","h","d","R","P0","Nb","a","Te","Pr"])    
    for i in range(len(BConfigMatrix)):
        for j in range(len(BConfigMatrix[i])):
            fileTreino.writerow(BConfigMatrix[i][j])
    print(".")

###########################################################################
##montar graficos
###########################################################################
            
#Gerar MATLAB
def gerarConjuntoMATLABR(BConfigMatrix):
    
    X = list(range(len(BConfigMatrix[0])))
    
    for j in range(len(BConfigMatrix[0])):
        X[j] = BConfigMatrix[0][j][3]
    
    for i in range(len(BConfigMatrix)):
        file = csv.writer(open("./matlab{0}a{1}h-R.csv".format(BConfigMatrix[i][0][6],BConfigMatrix[i][0][1]), "w"))
                
        for j in range(len(X)):
            file.writerow([BConfigMatrix[i][j][0],BConfigMatrix[i][j][3]])                        

#Gerar MATLAB
def gerarConjuntoMATLAB3D(BConfigMatrix):
    
    X = list(range(len(BConfigMatrix[0])))
    
    for j in range(len(BConfigMatrix[0])):
        X[j] = BConfigMatrix[0][j][3]
    
    for i in range(len(BConfigMatrix)):
        file = csv.writer(open("./matlab3D{0}a{1}h.csv".format(BConfigMatrix[i][0][6],BConfigMatrix[i][0][1]), "w"))
                
        for j in range(len(X)):
            file.writerow([BConfigMatrix[i][j][0],BConfigMatrix[i][j][3],BConfigMatrix[i][j][5]])    

def gerarConjuntoMATLABNb(BConfigMatrix):
    
    X = list(range(len(BConfigMatrix[0])))
    
    for j in range(len(BConfigMatrix[0])):
        X[j] = BConfigMatrix[0][j][5]
    
    for i in range(len(BConfigMatrix)):
        file = csv.writer(open("./matlab{0}a{1}h-Nb.csv".format(BConfigMatrix[i][0][6],BConfigMatrix[i][0][1]), "w"))
                
        for j in range(len(X)):
            file.writerow([BConfigMatrix[i][j][0],BConfigMatrix[i][j][5]]) 
            
##P0################
def montarGraphP0R(BConfigMatrix, Preal):
    legendaP0 = list(range(len(BConfigMatrix)))
    legendaP2 = list(range(len(BConfigMatrix)))
        
    P0  = list(range(len(BConfigMatrix)))
    P2 = list(range(len(BConfigMatrix)))
    
    R = list(range(len(BConfigMatrix[0])))
             
    #Formar Eihop
    i=0
    for i in range(len(BConfigMatrix)):
        legendaP0[i] = '{0} Salto(s) DNN'.format(BConfigMatrix[i][1][1])
        legendaP2[i] = '{0} Salto(s) MCE'.format(BConfigMatrix[i][1][1])
        
        #GERAR GRAFICOS PARA R
        if(len(BConfigMatrix[i])>1):
            auxP0 = list(range(len(BConfigMatrix[i])))
            auxP2 = list(range(len(BConfigMatrix[i])))
            
                                    
            for j in range(len(BConfigMatrix[i])):
                #auxP0[j] = math.log(BConfigMatrix[i][j][4],10)*10
                #auxP2[j] = math.log(Preal[i][j],10)*10
                auxP0[j] = BConfigMatrix[i][j][4]
                auxP2[j] = Preal[i][j]
                                            
            P0[i] = auxP0
            P2[i] = auxP2
            
    for j in range(len(BConfigMatrix[0])):
        R[j] = BConfigMatrix[0][j][3]
    
    legenda = []
    legenda.extend(legendaP0)
    legenda.extend(legendaP2)
    
    plotGraphP0(legenda, P0, P2, R, BConfigMatrix[0][0][5], True)

def montarGraphP0Nb(BConfigMatrix, Preal):
    legendaP0 = list(range(len(BConfigMatrix)))
    legendaP2 = list(range(len(BConfigMatrix)))
        
    P0  = list(range(len(BConfigMatrix)))
    P2 = list(range(len(BConfigMatrix)))
    
    Nb = list(range(len(BConfigMatrix[0])))
             
    #Formar Eihop
    i=0
    for i in range(len(BConfigMatrix)):
        legendaP0[i] = '{0} Salto(s) DNN'.format(BConfigMatrix[i][1][1])
        legendaP2[i] = '{0} Salto(s) MCE'.format(BConfigMatrix[i][1][1])
        
        #GERAR GRAFICOS PARA R
        if(len(BConfigMatrix[i])>1):
            auxP0 = list(range(len(BConfigMatrix[i])))
            auxP2 = list(range(len(BConfigMatrix[i])))
            
                                    
            for j in range(len(BConfigMatrix[i])):
                auxP0[j] = BConfigMatrix[i][j][4]
                auxP2[j] = Preal[i][j]
                                            
            P0[i] = auxP0
            P2[i] = auxP2
            
    for j in range(len(BConfigMatrix[0])):
        Nb[j] = BConfigMatrix[0][j][5]
    
    legenda = []
    legenda.extend(legendaP0)
    legenda.extend(legendaP2)
    
    plotGraphP0(legenda, P0, P2, Nb, BConfigMatrix[0][0][5], False)

#Funcao Plotar Graficos    
def plotGraphP0(legenda, P0, P2, X, a, isR):

    #MONTAR BRUNA
    plt.setp(plt.plot(X, P0[0], 'k3'), linewidth=(2))
    plt.setp(plt.plot(X, P0[1], 'cs'), linewidth=(1))
    plt.setp(plt.plot(X, P0[2], 'b*'), linewidth=(1))
    plt.setp(plt.plot(X, P0[3], 'rp'), linewidth=(1))
    
    #MONTAR PAULO
    plt.setp(plt.plot(X, P2[0], 'y--'), linewidth=(2))
    plt.setp(plt.plot(X, P2[1], 'g:'), linewidth=(3))
    plt.setp(plt.plot(X, P2[2], 'k-.'), linewidth=(2))
    plt.setp(plt.plot(X, P2[3], 'b-'), linewidth=(3))
    
    #Gerar figura
    plt.legend(legenda,loc='upper right')
    plt.ylabel('P0')
    
    if isR:
        plt.xlabel('R (bps)')
        plt.savefig('./R{0}-a{1}P0.eps'.format(X[len(X)-1],a), format='eps', dpi=1000)
    else:
        plt.xlabel('Nb (b)')
        plt.savefig('./Nb{0}-a{1}P0.eps'.format(X[len(X)-1],a), format='eps', dpi=1000)
    
    plt.show()        
  
def montarGraph2D(BConfigMatrix,x=3,y=0):
    
    legenda = list(range(len(BConfigMatrix)))
    Y = list(range(len(BConfigMatrix)))
    X = list(range(len(BConfigMatrix[0])))
             
    #Formar Eihop
    i=0
    for i in range(len(BConfigMatrix)):
        legenda[i] = '{0} Salto(s)'.format(BConfigMatrix[i][1][1])
        
        #GERAR GRAFICOS PARA R
        if(len(BConfigMatrix[i])>1):
            auxY = list(range(len(BConfigMatrix[i])))
                                    
            for j in range(len(BConfigMatrix[i])):
                auxY[j] = BConfigMatrix[i][j][y]
                                            
            Y[i] = auxY
    
    for j in range(len(BConfigMatrix[0])):
        X[j] = BConfigMatrix[0][j][x]
    
    plotGraph2D(legenda, Y, X, BConfigMatrix[0][0][5], x,y)

#Funcao Plotar Graficos    
def plotGraph2D(legenda, Y, X, a, x, y):
    
    lines = range(len(legenda))
        
    #EM FUNCAO DE R ou Nb
    for i in lines:
        plt.setp(plt.plot(X, Y[i]), linewidth=(i+1))
    
    #Gerar figura
    plt.legend(legenda,loc='upper right')
    plt.ylabel(definirLegenda(y))
    plt.xlabel(definirLegenda(x))

    plt.savefig('./2DGraph.a{0}.eps'.format(a), format='eps', dpi=1000)
         
    plt.show()

def montarGraph3D(BConfigMatrix,x=3,y=5,z=0):
    
    Z = list(range(len(BConfigMatrix)))
    X = list(range(len(BConfigMatrix)))
    Y = list(range(len(BConfigMatrix)))
    legenda = list(range(len(BConfigMatrix)))    
         
    #Formar Eihop3d
    for i in range(len(BConfigMatrix)):
        
        auxX = list(range(len(BConfigMatrix[i])))
        auxY = list(range(len(BConfigMatrix[i])))
        auxZ = list(range(len(BConfigMatrix[i])))        
        
        for j in range(len(BConfigMatrix[0])):
            auxX[j] = BConfigMatrix[i][j][x]#3R
            auxY[j] = BConfigMatrix[i][j][y]#5nb
            auxZ[j] = BConfigMatrix[i][j][z]#0EIHOP
        
        X[i] = auxX 
        Y[i] = auxY
        Z[i] = auxZ
        legenda[i] = '{0} Salto(s)'.format(BConfigMatrix[i][1][1])
        
    plotGraphEhop3D(legenda, Z, X, Y, BConfigMatrix[0][0][6],x,y,z)

#Funcao Plotar Graficos 3d    
def plotGraphEhop3D(legenda, Z, X, Y, a, x, y, z):
   
    lines = range(len(legenda))
    ax = plt.axes(projection="3d")
    
    for i in lines:
        plt.setp(ax.plot3D(X[i], Y[i], Z[i]))
    
    ax.set_xlabel(definirLegenda(x))
    ax.set_ylabel(definirLegenda(y))
    ax.set_zlabel(definirLegenda(z))
    
    #Gerar figura
    plt.legend(legenda,loc='upper right')    
    plt.savefig('./3DGraph.a{0}.eps'.format(a), format='eps', dpi=1000)

    plt.show()

#BEnh,h,d,R,P0,Nb,a,Te,Pr
def definirLegenda(a):
    switcher = {
        0: "Eihop (dBmJ/bit)",
        1: "Hops",
        2: "d (m)",
        3: "R (bps)",
        4: "Po (dBmJ)",
        5: "Nb (bits)",
        6: "alpha",
        7: "teta(0)"
    }
    return switcher.get(a)

#ENGENHARIA REVERSA BASE:
def baseTOmatrix():
    #CARREGAR BASE
    dataName = "/{0}.csv".format("CSMA.d80.NbO.PO")
    
    #ORDENA VALORES POR UMA DETERMINADA VARIÁVEL
    #orderBy = "h"
    #data = data.sort_values(by=[orderBy], ascending=True)
    
    print("Teste")

###########################################################################
##montar graficos VS
###########################################################################
def montarGraphDNNvs(BConfigMatrixDNN, BConfigMatrixMCE, _R, refGraph=0):
 
    #Formar Eihop - BRUNA
    EhopMCE = list(range(len(BConfigMatrixMCE)))
    EhopDNN = list(range(len(BConfigMatrixDNN)))
    for i in range(len(BConfigMatrixMCE)):
    
        #GERAR GRAFICOS PARA R
        if(len(BConfigMatrixMCE[i])>1):
            auxMCE = list(range(len(BConfigMatrixMCE[i])))
            auxDNN = list(range(len(BConfigMatrixDNN[i])))
                                    
            for j in range(len(BConfigMatrixMCE[i])):
                auxMCE[j] = BConfigMatrixMCE[i][j][refGraph]
                auxDNN[j] = BConfigMatrixDNN[i][j][refGraph]
                                            
            EhopMCE[i] = auxMCE
            EhopDNN[i] = auxDNN
    
    legendaDNN = ["1 hop(s) NN", "2 hop(s) NN", "4 hop(s) NN", "8 hop(s) NN"]
    legendaMCE = ["1 hop(s) ECM", "2 hop(s) ECM", "4 hop(s) ECM", "8 hop(s) ECM"]
    	
    plotGraphDNNvs(legendaMCE, EhopMCE, legendaDNN, EhopDNN, _R,refGraph)
    

def plotGraphDNNvs(legendaMCE, EhopMCE, legendaDNN, EhopDNN, _R, refGraph):
        
    x = _R
    fonte = 20
    linewidth = 3
    loclegend = 'upper right'
    typeG = 'R(bps)'
    
    #LABEL
    plt.ylabel('Eihop(dBmJ/bit)', fontsize=fonte)
    
    #label P0
    if(refGraph==4):plt.ylabel('P0(W)', fontsize=fonte)
    
    plt.xlabel(typeG, fontsize=fonte)
    
    #MONTAR PAULO
    plt.setp(plt.plot(x, EhopDNN[0], 'y--'), linewidth=(linewidth+2))
    plt.setp(plt.plot(x, EhopDNN[1], 'g:'), linewidth=(linewidth +3))
    plt.setp(plt.plot(x, EhopDNN[2], 'k-.'), linewidth=(linewidth+2))
    plt.setp(plt.plot(x, EhopDNN[3], 'b-'), linewidth=(linewidth +3))
    
    #MONTAR BRUNA
    plt.setp(plt.plot(x, EhopMCE[0], 'k3'), linewidth=(linewidth +2))
    plt.setp(plt.plot(x, EhopMCE[1], 'cs'), linewidth=(linewidth +1))
    plt.setp(plt.plot(x, EhopMCE[2], 'b*'), linewidth=(linewidth +1))
    plt.setp(plt.plot(x, EhopMCE[3], 'rp'), linewidth=(linewidth +1))
          
    #LEGENDAS        
    legenda = []
    legenda.extend(legendaDNN)
    legenda.extend(legendaMCE)
    plt.legend(legenda, ncol=2, loc=loclegend, prop={'size': 12})#upper - lower
    
    #PLOT
    
    plt.savefig('./GRAFICOFINAL.png', format='png', dpi=1000)
    plt.show()
    
    
    
###RUM
if __name__ == '__main__':
    baseTOmatrix()
