'''
ESSA CLASSE É USADA PARA MANIPULACAO DE DADOS
'''

import pandas as pd
import csv
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
import seaborn as sb

#CRIAR DADOS COM CLASSE < MINCLASSE
def DATA_addMCE(data):        
    
    data = data.values
    data2 = []
        
    fileTreino = csv.writer(open("./{0}.csv".format(nameData[MCE]), "w"))
    fileTreino.writerow(["BEnh","h","d","R","P0","Nb","a","Pr"])  
    
    for i in range(len(data)):
        data2.append(list(data[i]))
        data2[i].extend([MCE])
        fileTreino.writerow(data2[i])
            
    print("##############SALVO")   

#ORDENA DADOS E SALVA EM OUTRA BASE
def DATA_order(data):        
    
    #ORDENA VALORES POR UMA DETERMINADA VARIÁVEL
    orderBy = "BEnh"
    data = data.sort_values(by=[orderBy], ascending=True)
    data = data.values
    data2 = []
        
    fileTreino = csv.writer(open("./{0}order.csv".format(nameData[MCE]), "w"))
    fileTreino.writerow(["BEnh","h","d","R","P0","Nb","a","Pr"])  
    
    for i in range(len(data)):
        data2.append(list(data[i]))
        fileTreino.writerow(data2[i])
            
    print("##############SALVO")

#CRIA E SALVA UM CONJUNTO DE AMOSTRA DE DADO A PARTIR DO PANDAS
def DATA_cut(data):
    
    #treino e teste 20p
    data_train, data_test = train_test_split(data, test_size=(int(cutPercent*len(data))))
        
    print("base size: ",data.size)
    print("cut  size: ",data_test.size)
    
    #VERIFICAR DADOS
    DATA_preview(data)
    
    data = data_test.values
    data2 = []
    
    fileTreino = csv.writer(open("./{0}cut.csv".format(nameData[MCE]), "w"))
    fileTreino.writerow(["BEnh","h","d","R","P0","Nb","a","Pr"])  
    
    for i in range(len(data)):
        data2.append(list(data[i]))
        fileTreino.writerow(data2[i])
            
    print("##############DATAs SALVO")

    
#GERA IMAGENS DE PREVIEW DA BASE 
def DATA_preview(data):
    
    #treino e teste 20p
    data_train, data_test = train_test_split(data, test_size=(int(cutPercent*len(data))))

    print("base size: ",data.size)
    print("cut  size: ",data_test.size)
    
    data.hist(bins=50, figsize=(20,15))
    plt.savefig('./DATA.eps', format='eps', dpi=1000)
    plt.show()

    data_train.hist(bins=50, figsize=(20,15))
    plt.savefig('./DATAtreino.eps', format='eps', dpi=1000)
    plt.show()
            
    data_test.hist(bins=50, figsize=(20,15))
    plt.savefig('./DATAteste.eps', format='eps', dpi=1000)
    plt.show()

    #Avalia a representatividade dos dados para as CLASSES
    plt.figure(figsize=[10,4])
    sb.heatmap(data.corr())
    plt.savefig('./DATAclass.eps', format='eps', dpi=1000)
    plt.show()
    
    print("##############IMGs SALVO")    
    
def DATA_normalization(data):
    #BEnh,h,d,R,P0,Nb,a,Pr    
    out = 'BEnh'
    Xdata = data.drop(labels = [out], axis = 1)
    
    #normPr = 0#(data[6] - 0)/(0 - 0)
    
    func = "normh  = (data[0] - "+str(min(Xdata['h'])) +")/data[0] +"+str(max(Xdata['h'])  - (min(Xdata['h'])) )+"\n"
    func+= "normd  = (data[1] - "+str(min(Xdata['d'])) +")/data[1]"+str(max(Xdata['d'])  - (min(Xdata['d'])) )+"\n"
    func+= "normR  = (data[2] - "+str(min(Xdata['R'])) +")/data[2]"+str(max(Xdata['R'])  - (min(Xdata['R'])) )+"\n"
    func+= "normP0 = (data[3] - "+str(min(Xdata['P0']))+")/data[3]"+str(max(Xdata['P0']) - (min(Xdata['P0'])))+"\n"
    func+= "normNb = (data[4] - "+str(min(Xdata['Nb']))+")/data[4]"+str(max(Xdata['Nb']) - (min(Xdata['Nb'])))+"\n"
    func+= "norma  = (data[5] - "+str(min(Xdata['a'])) +")/data[5]"+str(max(Xdata['a'])  - (min(Xdata['a'])) )+"\n"
    func+= "normTe = 0 \n"
    #func+= "normTe = (data[6] - "+str(min(Xdata['Te']))+")/"+str(max(Xdata['Te']) - (min(Xdata['Te'])))+"\n" 
    func+= "data   = [normh,normd,normR,normP0,normNb,norma,normTe]\n"
    
    print(func)
    print((max(Xdata['h'])  - (min(Xdata['h']))))

def DATA_desnormalization(data):
    #BEnh,h,d,R,P0,Nb,a,Pr    
    out = 'BEnh'
    Xdata = data.drop(labels = [out], axis = 1)
    
    #normPr = 0#(data[6] - 0)/(0 - 0)
    
    func = "datah  = (norm[0]*("+str(max(Xdata['h']) -(min(Xdata['h'])))+")) + "+str(min(Xdata['h'])) +"\n"
    func+= "datad  = (norm[1]*("+str(max(Xdata['d']) -(min(Xdata['d'])))+")) + "+str(min(Xdata['d'])) +"\n"
    func+= "dataR  = (norm[2]*("+str(max(Xdata['R']) -(min(Xdata['R'])))+")) + "+str(min(Xdata['R'])) +"\n"
    func+= "dataP0 = (norm[3]*("+str(max(Xdata['P0'])-(min(Xdata['P0'])))+"))+ "+str(min(Xdata['P0']))+"\n"
    func+= "dataNb = (norm[4]*("+str(max(Xdata['Nb'])-(min(Xdata['Nb'])))+"))+ "+str(min(Xdata['Nb']))+"\n"
    func+= "dataa  = (norm[5]*("+str(max(Xdata['a']) -(min(Xdata['a'])))+")) + "+str(min(Xdata['a'])) +"\n"
    func+= "dataPr = 0 \n"
    #func+= "normTe = (data[6] - "+str(min(Xdata['Te']))+")/"+str(max(Xdata['Te']) - (min(Xdata['Te'])))+"\n" 
    func+= "data = [datah,datad,dataR,dataP0,dataNb,dataa,dataPr]\n"
    
    print(func)
    print((max(Xdata['h'])  - (min(Xdata['h']))))

#ORDENA DADOS E SALVA EM OUTRA BASE
def Pandas_view(data):        
    
    #ORDENA VALORES POR UMA DETERMINADA VARIÁVEL
    data = data.query('a == 4.0')
    #data = data.query('h == 1.0')
    data = data.query('d == 80')
    #data = data.query('BEnh < -150')    
        

    print(data)
        
    #print("##############SALVO")

##################################################################
###RUM
nameData = ["ALOHA", "CSMA"]
MCE = 1 #0 = ALOHA; 1 = CSMA; 2 = ALOHARM; 3 = CSMARM
dataName = "/{0}.csv".format(nameData[MCE])

MCE1 = 1
dataName1 = "/{0}.csv".format(nameData[MCE1])

cutPercent = 0.2
    
if __name__ == '__main__':
    
    #DATA    
    data0 = pd.read_csv("."+dataName)
    #data1 = pd.read_csv("."+dataName1)
    
    #concatenar bases
    #frames = [data0,data1]
    #data0 = pd.concat(frames)
    
    ###MANIPULACAO DE DADOS          
    #DATA_addMCE(data0) #CRIAR DADOS COM CLASSE < MINCLASSE
    #DATA_order(data0)  #ORDENA DADOS E SALVA EM OUTRA BASE
    #DATA_cut(data0)  #CRIA E SALVA UM CONJUNTO DE AMOSTRA DE DADO A PARTIR DO PANDAS  
    #DATA_preview(data0) #GERA IMAGENS DE PREVIEW DA BASE
    Pandas_view(data0)
    
    ###NORMALIZADOR MANUAL DE DADOS
    #DATA_normalization(data0)
    #DATA_desnormalization(data0)