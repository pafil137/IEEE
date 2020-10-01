'''
Um  Estudo  sobre  o  Consumo  de  Energia  em  RedesAd  HocLineares  Aloha  com  Saltos  Equidistantes
Funcao: Otimizar energia
ENTRADA:potencia  (pt),  tamanho  do  pacote(_Nb),  taxa (_R),  distancia (_dist),  numero  de  saltos (_h), Numero de Nos(h+1).
SAIDA: energia  gasta
'''

#biblioteca para Matematica Python
import math as math
import time

#UTil
import util_MCE  as util_MCE
#import ModeloEnergia.util_MCE  as util_MCE

#########ENTRADAS##############
####VARIAVEIS####
Pr = 0
nameProtocol = "ALOHA"

##Coeficiente  de  atenuacao
_a = [4]
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

_R2 = [1000]

#Tamanho  do  pacote (1bytes-1mb)
#_Nb = [1000]
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

#ppercentagem
#_Ppercent = 0.9
#_Pcount = 0.0001

_Ppercent = 0.1
_Pcount = 0.1

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

def functionPowerR(_h =_h, _dist=_dist, _R=_R, _Nb=_Nb, _a=_a):
    
    #Guardar melhor Configuracao
    _BConfigMatrix = list(range(len(_dist)*len(_h)*len(_a))) 
    _Preal = list(range(len(_dist)*len(_h)*len(_a)))
    
    i=0
    j=0
    
    inicio = time.time()
    
    for a in _a:        
        for dist in _dist:
            print("####A={0}-Dist={1}:R####".format(a, dist))
                        
            #Numero de Hops
            for h in _h:        
                print("#h{0}".format(h))
                #energia  gasta
                _BEnh = math.pow(10,7)
                _BConfigR = list(range(len(_R)))
                _Ptemp = list(range(len(_R)))
                
                #Contador R
                j=0
                
                #Taxa de transmissao
                for R in _R:
                    #print("#R{0}".format(R))
                    #criar isR: true:
                    if(len(_R)>1):
                        _BEnh = math.pow(10,7)
                        
                    #Tamanho do Pacote
                    for Nb in _Nb:
                        
                        #Distancia/numero de saltos
                        n = h + 1
                        d = dist/h
                                
                        muD = math.pow(((n-1)/float(n)),(n-2))
                        #mu = n/muD
                        mu = muD/n
                        
                        c1 = am*math.pow(d,a)*R*(1+Nb)
                        c2 = 4*bm*C2
                        c3 = math.pow(bamp,2)*math.pow(am,2)*math.pow(d,(2*a))*math.pow((1+Nb),2)
                        c4 = 8*bm*C2*math.pow(R,-1)*((C3*R)+(C4*Nb))*bamp*am*math.pow(d,a)
                        c5 = 4*bm*C2*bamp
                        
                        P0timo = (c1/c2) + R*(math.sqrt(c3+c4)/c5)
                        
                        #MONTAR RANGE DE P0
                        #contador de P0
                        Pvariation = _Ppercent*P0timo
                        Pcount = _Pcount*P0timo
                        
                        P0o = P0timo + Pvariation
                        P0i = P0timo - Pvariation
                                                
                        if(P0i<0):
                            P0i = Pcount                                             
                        
                        P0 = P0i                        
                        while(P0<P0o):                                                      
                            c6 = 1-(am*R*(math.pow(d,a))/float(2*bm*C2*P0))
                            
                            numeradorEnh = C3*R + (C4+bamp*P0)*Nb
                            denominadorEnh = R*Nb*math.pow(c6,Nb)
                            
                            ##########SAIDAS#############
                            #Energia dissipada em transmissao N hop
                            Enh= (numeradorEnh/denominadorEnh)*mu*h #%J/bit
                            Enh_dbm = 10*math.log10(Enh/0.001) #%dBmJ/bit
                            
                            ####SALVAR DADOS####
                            if(_BEnh > Enh_dbm):                    
                                _BEnh = Enh_dbm 
                                _BConfigR[j] = [_BEnh,h,d,R,math.log(P0,10)*10,Nb,a, Pr]
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
            print("####A={0}-Dist={1}:Nb####".format(a, dist))
                        
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
                    #print("#R{0}".format(R))
                    #criar isR: true:
                    if(len(_Nb)>1):
                        _BEnh = math.pow(10,7)
                        
                    #Tamanho do Pacote
                    for R in _R:
                        
                        #Distancia/numero de saltos
                        n = h + 1
                        d = dist/h
                                
                        muD = math.pow(((n-1)/float(n)),(n-2))
                        #mu = n/muD
                        mu = muD/n
                        
                        c1 = am*math.pow(d,a)*R*(1+Nb)
                        c2 = 4*bm*C2
                        c3 = math.pow(bamp,2)*math.pow(am,2)*math.pow(d,(2*a))*math.pow((1+Nb),2)
                        c4 = 8*bm*C2*math.pow(R,-1)*((C3*R)+(C4*Nb))*bamp*am*math.pow(d,a)
                        c5 = 4*bm*C2*bamp
                        
                        P0timo = (c1/c2) + R*(math.sqrt(c3+c4)/c5)
                        
                        #MONTAR RANGE DE P0
                        #contador de P0
                        Pvariation = _Ppercent*P0timo
                        Pcount = _Pcount*P0timo
                        
                        P0o = P0timo + Pvariation
                        P0i = P0timo - Pvariation
                                                
                        if(P0i<0):
                            P0i = Pcount                                             
                        
                        P0 = P0i                        
                        while(P0<P0o):                                                      
                            c6 = 1-(am*R*(math.pow(d,a))/float(2*bm*C2*P0))
                            
                            numeradorEnh = C3*R + (C4+bamp*P0)*Nb
                            denominadorEnh = R*Nb*math.pow(c6,Nb)
                            
                            ##########SAIDAS#############
                            #Energia dissipada em transmissao N hop
                            Enh= (numeradorEnh/denominadorEnh)*mu*h #%J/bit
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

        
#******************************PHYTHON FUNCTIONS*******************************
#************FUNCAO MAIN
if __name__ == '__main__':
    print("#########PARAMETROS##########\nEnergia dissipada-dBmJ/bit(Eih_DBM)\n *numero  de  saltos(i)\n distancia(d)\n taxa(R)\n  potencia(P)\n tamanho  do  pacote(Nb)\n Coeficiente  de  atenuacao (a)\n#############################\n")
    print("[Eih_dbm.0, h.1, d.2, R.3, P0.4, Nb.5, a.6, teta.7, Pr]")
      
    BConfigMatrixR0, PrealR = functionPowerR(_h, _dist, _R, _Nb, _a)
    BConfigMatrixR1, PrealR = functionPowerR(_h, _dist, _R, [1000], _a)
    
    util_MCE.montarGraphDNNvs(BConfigMatrixR0, BConfigMatrixR1,_R,4)
    
    #BConfigMatrixNb, PrealNb = functionPowerNb()
    
    if(len(_a)<2 and len(_dist)<2):
        #BConfigMatrixNb, PrealNb = functionPowerNb()
        print("montar3DGraph")
        #util_MCE.montarGraph3D(BConfigMatrixR)#x,y,z
        
        print("montarGraph2D")
        #util_MCE.montarGraph2D(BConfigMatrixR)
        
        print("montarGraphP0-R")
        #util_MCE.montarGraphP0R(BConfigMatrixR, PrealR)
               
        print("montarGraph2D-NB")
        #util_MCE.montarGraph2D(BConfigMatrixNb,5)
        
        print("montarGraphP0-NB")
        #util_MCE.montarGraphP0Nb(BConfigMatrixNb, PrealNb)
    
        ##R- TAXA###################
        print("montarDATAMatlab-R")
        #util_MCE.gerarConjuntoMATLABR(BConfigMatrixR)
        
        ##Nb- tamano pacote###################
        print("montarGraphMatlab-NB")
        #util_MCE.gerarConjuntoMATLABNb(BConfigMatrixNb)
        
    #print("gerarConjuntoTreino")
    #nameBase =nameProtocol+"{0}R".format(_dist[len(_dist)-1])
    
    
    #util_MCE.gerarConjuntoTreino(BConfigMatrixR, nameBase)
    
#    nameBase =nameProtocol+"{0}Nb".format(_dist[len(_dist)-1])
#    util_MCE.gerarConjuntoTreino(BConfigMatrixNb, nameBase)
    
    print("####END####")