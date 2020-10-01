'''
Um  Estudo  sobre  o  Consumo  de  Energia  em  RedesAd  HocLineares  Aloha  com  Saltos  Equidistantes
Funcao: Otimizar energia
ENTRADA:potencia  (pt),  tamanho  do  pacote(_Nb),  taxa (_R),  distancia (_dist),  numero  de  saltos (_h), Numero de Nos(h+1).
SAIDA: energia  gasta
'''

#biblioteca para Matematica Python
import math as math
import time
from numba import njit, prange
import numpy as np

#UTil
import util_MCE

#########ENTRADAS##############
####VARIAVEIS####
Pr = 3 #ALOHA RANDOM CODE

##Coeficiente  de  atenuacao
#_a = [2,4]
_a = range(2,6)

#Distancia  entre  fonte  e  destino
#_dist = [80]
_dist = range(1, 251, 1)
#_dist = range(30, 250, 5)

#Taxa de Transmissao (1bps - 1Mbps): 
#*quanto mais melhor
R0 = np.arange(10,int(math.pow(10,3)),10)
R1 = np.arange(int(math.pow(10,3)),int(math.pow(10,4)),int(math.pow(10,1)))
R2 = np.arange(int(math.pow(10,4)),int(math.pow(10,5)),int(math.pow(10,2)))
R3 = np.arange(int(math.pow(10,5)),int(math.pow(10,6)),int(math.pow(10,3)))
R4 = np.arange(int(math.pow(10,6)),int(math.pow(10,7)),int(math.pow(10,4)))
R5 = np.arange(int(math.pow(10,7)),int(math.pow(10,8)),int(math.pow(10,5)))
_R = np.concatenate((R0,R1,R2,R3,R4,R5), axis=None)
#_R = np.arange(50, int(math.pow(10,6)), 10000)

#Tamanho  do  pacote (1bytes-1mb)
#_Nb = [1000]
#_Nb = [1, 10, math.pow(10,2), math.pow(10,3), math.pow(10,4), math.pow(10,5), math.pow(10,6)]
_Nb = np.arange(1, 7981)#1kbp-8mb

#Numero de Nodos
#_n = range(2,101)
_n = [2, 5, 20, 70]

#teta(angulo de abertura)[rad][5 graus]
pi = math.pi 
_teta = np.array([0.15*(pi),0.30*(pi),0.45*(pi),0.60*(pi)])
#_teta = [0.15*(pi)]

####CONSTANTES####
#Potencia  do  circuito  de  recepcao
PrxElec = 279*math.pow(10,-3)
#Potencia  para  inicializacao
Pstart = 58.7*math.pow(10,-3)
#Tempo  de  inicializacao
Tstart = 446*math.pow(10,-6)
#Potencia  do  circuito  de  transmissao
PtxElec = 151*math.pow(10,-3)
#Nivel  de  potencia  Eq.  (4)

aamp = 174*math.pow(10,-3)
#Constante  de  proporcionalidade  Eq.  (4)
bamp = 5
#Densidade  espetral  de  ruido (Watt)
N0 = math.pow(10, -15.4)*(math.pow(10,-3))
#Frequencia  da  portadora
fc = 2.4*math.pow(10,9)
#Velocidade  da  luz
c = 3*math.pow(10,8)
#Ganho  da  antena  de  transmissao
Gtant = 1
#Ganho  da  antena  de  recepcao
Grant = 1
#Constante de Modulacao(BPSK)
am = 1
bm = 2
#variaveis nao identificadas
L = 1

####FORMULAS#####
#comprimento de onda (lambda)
lamb = (c/fc)

###consumo  medio  de  energia  por  bit  para  umatransmissao  por  multiplos  saltos
##Ec1 (Formula 6.a - pg2)
C3 = (2*Tstart*Pstart) 

##Ec (Formula 6.b - pg2)
C4 = PtxElec+PrxElec+aamp

##c2 (Formula 10 - pg3)
C2 = (Gtant*Grant*math.pow(lamb,2))/(math.pow(4*math.pi,2)*N0*L)

def functionPower():
    
    #Guardar melhor Configuracao
    _BConfigMatrix = list(range(len(_dist)*len(_n)*len(_a))) 
    
    inicio = time.time()
    i = 0
    
    for a in _a:        
        for dist in _dist:
            print("####A={0}-Dist={1}####".format(a, dist))
                        
            #Numero de Nodos
            for n in _n:
                print("#n{0}".format(n))
                _BConfigMatrix[i] = GPU(a, dist, n)
                i+=1
    
    fim = time.time()
    print("tempo:", fim - inicio)
    
    return _BConfigMatrix

#CUDA PARALELO GPU
@njit(parallel=True)
def GPU(a,dist,n):

    nR = _R.shape[0]
    nNb = _Nb.shape[0]
    nTeta = _teta.shape[0]
    
    #tamanho vetor saida = 9
    s = (len(_R),9)
    _BConfigR = np.zeros(s)
    
    #Taxa de transmissao
    for i in prange(nR):
        
        R = _R[i]
        BEnh = math.pow(10,7)
        
        #Tamanho do Pacote
        for j in prange(nNb):
            Nb = _Nb[j]
            
            #angulo de busca do nodo
            for t in prange(nTeta):
                teta = _teta[t]
                
                #Distancia/numero de saltos
                #Numero  de  saltos  ate  o  destino
                h = int(((dist/3)*(math.pow((math.sqrt(2)+math.log(1+math.sqrt(2))),2)))/(math.sqrt((pi*math.pow(dist,2)/(2*n*teta)))*(2/teta)*(math.sin(teta/2))))            
                        
                d = (math.sqrt(pi/(2*(n/(math.pow(dist,2)))*teta)));
                                    
                muD = math.pow(((n-1)/float(n)),(n-2))
                
                mu = muD/n
                
                c1 = am*math.pow(d,a)*R*(1+Nb)
                c2 = 4*bm*C2
                c3 = math.pow(bamp,2)*math.pow(am,2)*math.pow(d,(2*a))*math.pow((1+Nb),2)
                c4 = 8*bm*C2*math.pow(R,-1)*((C3*R)+(C4*Nb))*bamp*am*math.pow(d,a)
                c5 = 4*bm*C2*bamp
                
                P0 = (c1/c2) + R*(math.sqrt(c3+c4)/c5)
                                                                      
                c6 = 1-(am*R*(math.pow(d,a))/float(2*bm*C2*P0))
                
                numeradorEnh = C3*R + (C4+bamp*P0)*Nb
                denominadorEnh = R*Nb*math.pow(c6,Nb)
                
                ##########SAIDAS#############
                #Energia dissipada em transmissao N hop
                Enh= (numeradorEnh/denominadorEnh)*mu*h #%J/bit
                Enh_dbm = 10*math.log10(Enh/0.001) #%dBmJ/bit
                
                ####SALVAR DADOS####
                if(BEnh > Enh_dbm):                    
                    BEnh = Enh_dbm 
                    _BConfigR[i] = [BEnh,h,d,R,P0,Nb,a, teta, Pr]
        
    return _BConfigR
        
#******************************PHYTHON FUNCTIONS*******************************
#************FUNCAO MAIN
if __name__ == '__main__':
    print("#########PARAMETROS##########\nEnergia dissipada-dBmJ/bit(Eih_DBM)\n *numero  de  saltos(i)\n distancia(d)\n taxa(R)\n  potencia(P)\n tamanho  do  pacote(Nb)\n Coeficiente  de  atenuacao (a)\n#############################\n")
    print("[Eih_dbm, i, d, R, P0, Nb, a, Pr]")
    
    BConfigMatrix = functionPower()
        
    print("gerarConjuntoTreino")
    nameBase ="ALOHA"+"{0}".format(_R[len(_R)-1])
    util_MCE.gerarConjuntoTreino(BConfigMatrix, nameBase)
    
    print("####END####")
