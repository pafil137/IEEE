'''
Um  Estudo  sobre  o  Consumo  de  Energia  em  RedesAd  HocLineares  Aloha  com  Saltos  Equidistantes
Funcao: Otimizar energia
ENTRADA:potencia  (pt),  tamanho  do  pacote(_Nb),  taxa (_R),  distancia (_dist),  numero  de  saltos (_h), Numero de Nos(h+1).
SAIDA: energia  gasta
'''

#biblioteca para Matematica Python
import math as math
import time
import sys

#UTil
import util_MCE

#import _pickle as cPickle

#########ENTRADAS##############
####VARIAVEIS####
Pr = 1
nameProtocol = "CSMA"

##Coeficiente  de  atenuacao
_a = [2]
#_a = [2,4]

#Distancia  entre  fonte  e  destino
_dist = [80]
#_dist = [30,80,150,250]
#_dist = range(10, 250, 5)

#Taxa de Transmissao (1bps - 1Mbps): 
#*quanto mais melhor

_R = []
_R.extend(range(int(math.pow(10,3)),int(math.pow(10,4)),int(math.pow(10,3))))
_R.extend(range(int(math.pow(10,4)),int(math.pow(10,5)),int(math.pow(10,4))))
_R.extend(range(int(math.pow(10,5)),int(math.pow(10,6)),int(math.pow(10,5))))
_R.extend(range(int(math.pow(10,6)),int(math.pow(10,7)),int(math.pow(10,6)))) 

#_R = range(50, int(math.pow(10,6)), 10000)
#_R = range(50, int(math.pow(10,6)), 50)

_R2 = [1000000]

#Tamanho  do  pacote (1bytes-1mb)
#_Nb = [1000]
#_Nb = [1, 10, math.pow(10,2), math.pow(10,3), math.pow(10,4), math.pow(10,5), math.pow(10,6)]
_Nb = range(1, 7981)#1kbp-8mb

_Nb2 = []
_Nb2.extend(range(int(math.pow(10,1)),int(math.pow(10,2)),int(math.pow(10,1))))
_Nb2.extend(range(int(math.pow(10,2)),int(math.pow(10,3)),int(math.pow(10,2))))
_Nb2.extend(range(int(math.pow(10,3)),int(math.pow(10,4)),int(math.pow(10,3))))
_Nb2.extend(range(int(math.pow(10,4)),int(math.pow(10,5)),int(math.pow(10,4))))


#Numero  de  saltos  ate  o  destino
#_h = [1]
_h = [1, 2, 4, 8]

#Numero de Nodos
n=0
#_n = h + 1

#ppercentagem
_Pi = 0.01#0.001
_Pf = 0.1
Pcount = _Pi

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

def functionPowerR(_h =_h, _dist=_dist, _R=_R, _Nb=_Nb, _a=_a, p0type='Pa'):
    
    #Guardar melhor Configuracao
    _BConfigMatrix = list(range(len(_dist)*len(_h)*len(_a))) 
    _Preal = list(range(len(_dist)*len(_h)*len(_a)))
    
    i=0
    j=0
    
    inicio = time.time()
    
    for a in _a:
                
        for dist in _dist:
            print("####A={0}-Dist={1}####".format(a, dist))
                        
            #Numero de Hops
            for h in _h:        
                ###(OK)### numero de saltos
                print("#h{0}".format(h))
                
                ###(OK)### tau
                n = h + 1
                tau = tauFind(n)
                    
                print("h:",h)
                
                print("#n{0}".format(h))
                #energia  gasta
                _BEnh = math.pow(10,7)
                _BConfigR = list(range(len(_R)))
                _Ptemp = list(range(len(_R)))
                
                #Contador R
                j=0
                
                #Taxa de transmissao
                for R in _R:                    
                    #criar isR: true:
                    if(len(_R)>1):
                        _BEnh = math.pow(10,7)
                        
                    #Tamanho do Pacote
                    for Nb in _Nb:
                                            
                        #Distancia/numero de saltos                        
                        d = dist/h
                                
                        ###(OK)### Bianchi(10)
                        Ptr = 1 - math.pow((1 - tau),n)
                        
                        ###(OK)### Bianchi(11) probabilidade de sucesso csma/ca
                        ps  = (n*tau*math.pow((1-tau),(n - 1)))/Ptr
                        
                        ###(OK)###                        
                        c1 = am*math.pow(d,a)*R*(1+Nb) #(OK)
                        c2 = 4*bm*C2 #(OK)
                        c3 = math.pow(bamp,2)*math.pow(am,2)*math.pow(d,(2*a))*math.pow((1+Nb),2)#(OK)
                        c4 = 8*bm*C2*math.pow(R,-1)*((C3*R)+C4*Nb)*bamp*am*math.pow(d,a)#(OK)
                        c5 = 4*bm*C2*bamp#(OK)
                        
                        P0timo = (c1/c2) + R*math.sqrt(c3+c4)/c5#(OK)
                        
                        ###MONTAR P0
                        #BRUNA       
                        #BALIZADO
                        P0i = P0timo - 0.9*P0timo
                        P0f = P0timo + 0.9*P0timo
                        
                        if (p0type=='Br'):
                            P0i = P0timo#_Pi
                            P0f = P0timo#_Pf
                        
                        #PAULO
                        #P0i = _Pi 
                        #P0f = _Pf                                                                                             
                        
                        P0 = P0i                                                                       
                        while(P0<=P0f):   
                            Prs = math.pow(1- ((am*math.pow(d,a)*R)/(float(2*bm*C2*P0))),Nb) #(OK)                             
                            numeradorEnh = C3*R + (C4 + bamp*P0)*Nb #(OK) 
                            denominadorEnh = R*Nb*Prs #(OK)
                            
                            ##########SAIDAS#############                            
                            Enh= (numeradorEnh/denominadorEnh)*h*(1/ps) #%J/bit
                            Enh_dbm = 10*math.log10(Enh/0.001) #%dBmJ/bit
                                                                
                            ####SALVAR DADOS####
                            if(_BEnh > Enh_dbm):                    
                                _BEnh = Enh_dbm 
                                _BConfigR[j] = [_BEnh,h,d,R,P0,Nb,a, 0, Pr]
                                _Ptemp[j] = P0timo
                                
                            P0+=Pcount
                             
                    j+=1
                    
                _BConfigMatrix[i] = _BConfigR
                _Preal[i] = _Ptemp
                i+=1
    
    fim = time.time()
    print("tempo:", fim - inicio)
    
    return _BConfigMatrix, _Preal

def functionPowerNb(_h =_h, _dist=_dist, _R=_R2, _Nb=_Nb2, _a=_a):
    
    #Guardar melhor Configuracao
    _BConfigMatrix = list(range(len(_dist)*len(_h)*len(_a)))
    _Preal = list(range(len(_dist)*len(_h)*len(_a)))
    
    i=0
    j=0
    
    inicio = time.time()
    
    for a in _a:
                
        for dist in _dist:
            print("####A={0}-Dist={1}####".format(a, dist))
                        
            #Numero de Hops
            for h in _h:        
                print("#h{0}".format(h))
                #energia  gasta
                _BEnh = math.pow(10,7)
                _BConfigR = list(range(len(_Nb)))
                _Ptemp = list(range(len(_Nb)))
                
                #Contador R
                j=0
                
                #Taxa de transmissao
                for Nb in _Nb:
                    
                    #criar isR: true:
                    if(len(_Nb)>1):
                        _BEnh = math.pow(10,7)
                        
                    #Tamanho do Pacote
                    for R in _R:
                                            
                        #Distancia/numero de saltos
                        n = h + 1
                        d = dist/h
                                
                        tau = tauFind(n)
                        Ptr     = 1 - math.pow((1 - tau),n); #Bianchi(10)
                        ps      = (n*tau*math.pow((1-tau),(n - 1)))/Ptr; #Bianchi(11) - probabilidade de sucesso csma/ca
                                                
                        c1 = am*math.pow(d,a)*R*(1+Nb)
                        c2 = 4*bm*C2
                        c3 = math.pow(bamp,2)*math.pow(am,2)*math.pow(d,(2*a))*math.pow((1+Nb),2)
                        c4 = 8*bm*C2*math.pow(R,-1)*((C3*R)+(C4*Nb))*bamp*am*math.pow(d,a)
                        c5 = 4*bm*C2*bamp
                        
                        P0timo = (c1/c2) + R*(math.sqrt(c3+c4)/c5)
                        
                        ###MONTAR P0
                        #BRUNA       
                        #P0i = P0timo#_Pi
                        #P0f = P0timo#_Pf
                        
                        #BALIZADO
                        P0i = P0timo - 0.9*P0timo
                        P0f = P0timo + 0.9*P0timo
                        
                        #PAULO
                        #P0i = _Pi 
                        #P0f = _Pf                                                                                             
                        
                        P0 = P0i                          
                        while(P0<P0f):   
                                                                                        
                            c6 = 1-(am*R*(math.pow(d,a))/float(2*bm*C2*P0))    
                            
                            numeradorEnh = C3*R + (C4+bamp*P0)*Nb
                            denominadorEnh = R*Nb*math.pow(c6,Nb)
                            
                            ##########SAIDAS#############
                            
                            #inconcistencia matematica quando R=2500/(Enh=negativo)
                            Enh= (numeradorEnh/denominadorEnh)*(1/float(ps))*h #%J/bit
                            prs = ps*math.pow(c6,Nb)
                            EC = 0
                            
                            for k in range(7):
                                EC += k*math.pow(1- prs,k-1)
                            
                            #Energia dissipada em transmissao N hop    
                            E1hop = prs*numeradorEnh*EC/(R*Nb) 
                            Enh= E1hop*h #%J/bit
                            Enh_dbm = sys.float_info.max
                                     
                            if(Enh > 0):
                                Enh_dbm = 10*math.log10(Enh/0.001) #%dBmJ/bit
                                                                
                            ####SALVAR DADOS####
                            if(_BEnh > Enh_dbm):                    
                                _BEnh = Enh_dbm 
                                _BConfigR[j] = [_BEnh,h,d,R,P0,Nb,a, Pr]
                                _Ptemp[j] = P0timo
                                
                            P0+=Pcount    
                    j+=1
                    
                _BConfigMatrix[i] = _BConfigR
                _Preal[i] = _Ptemp
                i+=1
    
    fim = time.time()
    print("tempo:", fim - inicio)
    
    return _BConfigMatrix, _Preal

#Funcao encontrar tal
def tauFind(n):
    Ret = 7
    
    #p=0.1
    thistau = 0.1
    
    p = math.pow((1 - (1-thistau)),(n-1))
    thistau = 1/(1 + ((1-p)/(2*(1-math.pow(p,(Ret+1)))))*(findSomatorio(Ret,p)))
    print(thistau)
    return thistau

#%Funcao do tamanho da janela de backoff (2^5) - throughput e praticamente constante; bianchi1 2000
def findSomatorio(Ret, p):
    somatorio = 0
    W = 32  
    for j in range(Ret):
        somatorio = somatorio + (math.pow(p,j))*(math.pow(2,j)*W-1) - (1-math.pow(p,(Ret+1)))
    
    return somatorio
        
#******************************PHYTHON FUNCTIONS*******************************
#************FUNCAO MAIN
if __name__ == '__main__':
    print("#########PARAMETROS##########\nEnergia dissipada-dBmJ/bit(Eih_DBM)\n *numero  de  saltos(i)\n distancia(d)\n taxa(R)\n  potencia(P)\n tamanho  do  pacote(Nb)\n Coeficiente  de  atenuacao (a)\n#############################\n")
    print("[Eih_dbm.0, h.1, d.2, R.3, P0.4, Nb.5, a.6, teta.7, Pr]")
      
    #BConfigMatrixR, PrealR = functionPowerR()
    #BConfigMatrixNb, PrealNb = functionPowerNb()
    
    BConfigMatrixR0, PrealR = functionPowerR()
    BConfigMatrixR1, PrealR = functionPowerR(_h, _dist, _R, [1000], _a, 'Br')
    #BConfigMatrixR0, PrealR = functionPowerR(_h, _dist, _R, [1000], _a, 'Br')
    
    
    util_MCE.montarGraphDNNvs(BConfigMatrixR0, BConfigMatrixR1,_R,4)
    
    if(len(_a)<2 and len(_dist)<2):
        #BConfigMatrixNb, PrealNb = functionPowerNb()
        print("montar3DGraph")
        #util_MCE.montarGraph3D(BConfigMatrixR,3,5,4)#x,y,z
        
        print("montarGraph2D")
        #util_MCE.montarGraph2D(BConfigMatrixR)
        
        print("montarGraphP0-R")
        #util_MCE.montarGraphP0R(BConfigMatrixR, PrealR)              
    
        ##R- TAXA###################
        print("montarDATAMatlab-R")
        #util_MCE.gerarConjuntoMATLABR(BConfigMatrixR)
        
    print("gerarConjuntoTreino")
    #nameBase =nameProtocol+"{0}R".format(_dist[len(_dist)-1])
    #util_MCE.gerarConjuntoTreino(BConfigMatrixR, nameBase)
    
    #nameBase =nameProtocol+"{0}Nb".format(_dist[len(_dist)-1])
    #util_MCE.gerarConjuntoTreino(BConfigMatrixNb, nameBase)
    
    print("####END####")
