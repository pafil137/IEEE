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
Pr = 1 #CSMA CODE

##Coeficiente  de  atenuacao
_a = [4]
#_a = range(2,6)

#Distancia  entre  fonte  e  destino
#_dist = [40]
_dist = range(10, 251, 10)
#_dist = range(30, 250, 5)

#Taxa de Transmissao (1bps - 1Mbps): 
#*quanto mais melhor
R0 = np.arange(10,int(math.pow(10,3)),10)
R1 = np.arange(int(math.pow(10,3)),int(math.pow(10,4)),int(math.pow(10,1)))
R2 = np.arange(int(math.pow(10,4)),int(math.pow(10,5)),int(math.pow(10,2)))
R3 = np.arange(int(math.pow(10,5)),int(math.pow(10,6)),int(math.pow(10,3)))
R4 = np.arange(int(math.pow(10,6)),int(math.pow(10,7)),int(math.pow(10,4)))
R5 = np.arange(int(math.pow(10,7)),int(math.pow(10,8)),int(math.pow(10,5)))
#_R = np.concatenate((R0,R1,R2,R3,R4,R5), axis=None)
_R = np.arange(50, int(math.pow(10,7)), 10000)

#Tamanho  do  pacote (1bytes-1mb)
#_Nb = [1000]
#_Nb = [1, 10, math.pow(10,2), math.pow(10,3), math.pow(10,4), math.pow(10,5), math.pow(10,6)]
_Nb = np.arange(1, 7981)#1kbp-8mb

#Numero  de  saltos  ate  o  destino
_h = [1]
#_h = range(1,9)
#_h = [1, 2, 4, 8]

#Numero de Nodos
#_n = h + 1
n=0

#ppercentagem
_Pcount = 0.01

#teta(angulo de abertura)
teta = 0

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
#PI
pi = math.pi 
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
    _BConfigMatrix = list(range(len(_dist)*len(_h)*len(_a))) 
    
    inicio = time.time()
    i = 0
    
    for a in _a:        
        for dist in _dist:
            print("####A={0}-Dist={1}:R####".format(a, dist))
                        
            #Numero de Hops
            for h in _h:        
                print("#h{0}".format(h))
                                          
                _BConfigMatrix[i] = GPU(a, dist, h)
                i+=1
    
    fim = time.time()
    print("tempo:", fim - inicio)
    
    return _BConfigMatrix

#CUDA PARALELO GPU
@njit(parallel=True)
def GPU(a,dist,h):

    nR = _R.shape[0]
    nNb = _Nb.shape[0]
    
    s = (len(_R),9)
    BConfigR = np.zeros(s)
    
    #Taxa de transmissao
    for i in prange(nR):
        
        R = _R[i]
        BEnh = math.pow(10,7)
        
        #Tamanho do Pacote
        for j in prange(nNb):
            
            Nb = _Nb[j]
                                
            #Distancia/numero de saltos
            n = h + 1
            d = dist/h
            
            ##############TAU FIND
            Ret = 7
    
            #p=0.1
            thistau = 0.1
   
            p = 1 - math.pow(1-thistau,(n-1))
            
            somatorio = 0
            W = 32
            
            for k in prange(Ret):
                somatorio = somatorio + (math.pow(p,k))*(math.pow(2,k)*W-1) - (1-math.pow(p,(Ret+1)))
            
            thistau = 1/(1 + ((1-p)/(2*(1-math.pow(p,(Ret+1)))))*(somatorio))
            tau = thistau
                        
            Ptr     = 1 - math.pow((1 - tau),n); #Bianchi(10)
            ps      = (n*tau*math.pow((1-tau),(n - 1)))/Ptr; #Bianchi(11) - probabilidade de sucesso csma/ca
                                    
            c1 = am*math.pow(d,a)*R*(1+Nb)
            c2 = 4*bm*C2
            c3 = math.pow(bamp,2)*math.pow(am,2)*math.pow(d,(2*a))*math.pow((1+Nb),2)
            c4 = 8*bm*C2*math.pow(R,-1)*((C3*R)+(C4*Nb))*bamp*am*math.pow(d,a)
            c5 = 4*bm*C2*bamp
            
            P0timo = (c1/c2) + R*(math.sqrt(c3+c4)/c5)
            
            #MONTAR RANGE DE P0
            #contador de P0
            P0i = P0timo - 0.9*P0timo
            P0f = P0timo - 0.7*P0timo            
                        
            #P0 = P0i            
            _P = np.arange(P0i,P0f, _Pcount)            
            nP0 = _P.shape[0]
                        
            for l in prange(nP0):
                
                P0 = _P[l]
                                            
                c6 = 1-(am*R*(math.pow(d,a))/float(2*bm*C2*P0))    
                
                numeradorEnh = C3*R + (C4+bamp*P0)*Nb
                #denominadorEnh = R*Nb*math.pow(c6,Nb)
                #Enh= (numeradorEnh/denominadorEnh)*(1/float(ps))*h #%J/bit
                #Enh_dbm = math.pow(10,7)
                
                ##########SAIDAS#############                
                #inconcistencia matematica quando R=2500/(Enh=negativo)                
                prs = ps*math.pow(c6,Nb)
                EC = 0
                
                for k in prange(7):
                    EC += k*math.pow(1- prs,k-1)
                
                #Energia dissipada em transmissao N hop    
                E1hop = prs*numeradorEnh*EC/(R*Nb)
                Enh= E1hop*h #%J/bit
                         
                Enh_dbm = 10*math.log10(Enh/0.001) #%dBmJ/bit
            
                ####SALVAR DADOS####
                if(BEnh > Enh_dbm):                    
                    BEnh = Enh_dbm 
                    BConfigR[i] = [BEnh,h,d,R,P0,Nb,a, teta,Pr]
                #P0+=_Pcount 
            
    return BConfigR
        
#******************************PHYTHON FUNCTIONS*******************************
#************FUNCAO MAIN
if __name__ == '__main__':
    print("#########PARAMETROS##########\nEnergia dissipada-dBmJ/bit(Eih_DBM)\n *numero  de  saltos(i)\n distancia(d)\n taxa(R)\n  potencia(P)\n tamanho  do  pacote(Nb)\n Coeficiente  de  atenuacao (a)\n#############################\n")
    print("[Eih_dbm, i, d, R, P0, Nb, a, Pr]")
    
    BConfigMatrix = functionPower()
        
    print("gerarConjuntoTreino")
    nameBase ="CSMA"+"{0}".format(_dist[len(_dist)-1])
    util_MCE.gerarConjuntoTreino(BConfigMatrix, nameBase)
    
    print("####END####")
