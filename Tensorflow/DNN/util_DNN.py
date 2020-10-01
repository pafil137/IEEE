import tensorflow as tf
from tensorflow.lite.python import lite

from sklearn.metrics import mean_squared_error

import matplotlib.pyplot as plt
import seaborn as sb


'''
for i in range(len(Blist)):
    for j in range(len(Blist[i])):        
            print("Eihop", Blist[i][j][0])
            print("i", Blist[i][j][1])
            print("d", Blist[i][j][2])
            print("R", Blist[i][j][3])
            print("P0", Blist[i][j][4])
            print("Nb", Blist[i][j][5])
            print("a", Blist[i][j][6])
'''

#DNN-CONVERT############################################################################  
####NORMALIZATION
def normalizationMinMax0(Xdata):
    
    normPr = Xdata['Pr']
    normalized = []
    
    #BEnh,h,d,R,P0,Nb,a,Pr
    for data in Xdata.values:
        
        normh  = (data[0] - 1.0)/7.0
        normd  = (data[1] - 1.25)/243.75
        normR  = (data[2] - 50.0)/990000.0
        normP0 = (data[3] - 1.7071578126875538e-06)/10.694997679940673
        normNb = (data[4] - 1.0)/7979.0
        norma  = (data[5] - 2.0)/2.0
        normPr = (data[6]) 
        
        normalized.append([normh,normd,normR,normP0,normNb,norma,normPr])
        
    return normalized

###KERAS TO PD
def freeze_session(session, keep_var_names=None, output_names=None, clear_devices=True):
    """
    Freezes the state of a session into a pruned computation graph.
    Creates a new computation graph where variable nodes are replaced by
    constants taking their current value in the session. The new graph will be
    pruned so subgraphs that are not necessary to compute the requested
    outputs are removed.
    @param session The TensorFlow session to be frozen.
    @param keep_var_names A list of variable names that should not be frozen,
                          or None to freeze all the variables in the graph.
    @param output_names Names of the relevant graph outputs.
    @param clear_devices Remove the device directives from the graph for better portability.
    @return The frozen graph definition.
    """
    from tensorflow.python.framework.graph_util import convert_variables_to_constants
    graph = session.graph
    with graph.as_default():
        freeze_var_names = list(set(v.op.name for v in tf.global_variables()).difference(keep_var_names or []))
        output_names = output_names or []
        output_names += [v.op.name for v in tf.global_variables()]
        input_graph_def = graph.as_graph_def()
        if clear_devices:
            for node in input_graph_def.node:
                node.device = ""
        frozen_graph = convert_variables_to_constants(session, input_graph_def,
                                                      output_names, freeze_var_names)
        return frozen_graph

#DNN-PLOT############################################################################
def debugDNN(yt = None, pred= None, data= None, dataT= None, mse=[], nameData= ""):
    
    msePred = mean_squared_error(yt, pred)
    print('M.S.E:', msePred)
    
    fonte = 14
    
    if((not data.empty)and(not dataT.empty)):
        
        data.hist(bins=50, figsize=(20,15))
        plt.savefig('./OUTPUT/{0}dataTrain.eps'.format(nameData), format='eps', dpi=1000)
        plt.show()
    
                
        dataT.hist(bins=50, figsize=(20,15))
        plt.savefig('./OUTPUT/{0}dataTest.eps'.format(nameData), format='eps', dpi=1000)
        plt.show()
    
        #Avalia a representatividade dos dados para as CLASSES
        plt.figure(figsize=[10,4])
        sb.heatmap(data.corr())
        plt.savefig('./OUTPUT/{0}train.eps'.format(nameData), format='eps', dpi=1000)
        plt.show()
         
    #PLOTAR GRAFICO PREDICAO    
    #MSE - EVOLUTION
    if(len(mse)> 1):
        plt.ylabel('MSE (dBmJ/bit)\u00B2', fontsize=fonte)
        plt.xlabel('NN training interactions', fontsize=fonte)
        plt.plot(mse)       
        plt.savefig('./OUTPUT/{0}MSE.eps'.format(nameData), format='eps', dpi=1000)     
        plt.show()
            
    #PREDICTED VS RIGHT
    plt.ylabel('Eihop(dBmJ/bit)', fontsize=fonte)
    plt.xlabel('NN training interactions', fontsize=fonte)
    
    plt.plot(yt, 'ro')
    plt.plot(pred, "b+")
    plt.legend(["NN training","Eihop CSMA/CA Model"], prop={'size': fonte})
    plt.savefig('./OUTPUT/{0}POINT.eps'.format(nameData), format='eps', dpi=1000)    
    plt.show()
    
    #PREDICTED VS RIGHT
    plt.ylabel('Eihop(dBmJ/bit)', fontsize=fonte)
    plt.xlabel('Test Set Samples',fontsize=fonte)
    
    plt.plot(sorted(yt), "y-")
    plt.plot(sorted(pred), "k:")
    #plt.legend(["Eihop CSMA/CA Model","DNN training"], prop={'size': fonte})
    plt.legend(["Eihop "+nameData+" Model","NN training"], prop={'size': fonte})
    plt.savefig('./OUTPUT/{0}PvR.eps'.format(nameData), format='eps', dpi=1000)    
    plt.show()

def litedebugDNN(yt = None, pred= None, mse = [], nameData= "", nameTest = ""):
    
    msePred = mean_squared_error(yt, pred)
    print('M.S.E:', msePred)
    
    fonte = 14
         
    #PLOTAR GRAFICO PREDICAO    
    #MSE - EVOLUTION
    if(len(mse)> 1):
        plt.ylabel('MSE (dBmJ/bit)\u00B2', fontsize=fonte)
        plt.xlabel('NN training interactions', fontsize=fonte)
        plt.plot(mse)       
        plt.savefig('./OUTPUT/{0}MSE.eps'.format(nameData), format='eps', dpi=1000)     
        plt.show()
            
        #PREDICTED VS RIGHT
    plt.ylabel('Eihop(dBmJ/bit)', fontsize=fonte)
    plt.xlabel('Test Set Samples',fontsize=fonte)
    
    plt.plot(sorted(yt), "y-")
    plt.plot(sorted(pred), "k:")
    #plt.legend(["Eihop CSMA/CA Model","DNN training"], prop={'size': fonte})
    plt.legend(["Eihop "+nameTest+" Model","NN training"], prop={'size': fonte})
    plt.savefig('./OUTPUT/{0}PvR.eps'.format(nameData), format='eps', dpi=1000)    
    plt.show()


def DNNdataGeneration(BConfigMatrix):
    
    Y = []
    X = []
    
    #GERAR ENTRADA ESPECIFICAS DNN  
    for i in range(len(BConfigMatrix)):
        for j in range(len(BConfigMatrix[i])):
            Y.append(BConfigMatrix[i][j][0])
            X.append([BConfigMatrix[i][j][1],BConfigMatrix[i][j][2],BConfigMatrix[i][j][3],BConfigMatrix[i][j][4],BConfigMatrix[i][j][5],BConfigMatrix[i][j][6]])
            
    return X, Y

#Gerar MATLAB
def montarGraphDNNvs(BConfigMatrix, preDNN, _R, typeG):
     
    #Formar Eihop - BRUNA
    EhopMCE = list(range(len(BConfigMatrix)))
    for i in range(len(BConfigMatrix)):

        
        #GERAR GRAFICOS PARA R
        if(len(BConfigMatrix[i])>1):
            auxE = list(range(len(BConfigMatrix[i])))
                                    
            for j in range(len(BConfigMatrix[i])):
                auxE[j] = BConfigMatrix[i][j][0]
                                            
            EhopMCE[i] = auxE
    
    #Formar Eihop - PAULO
    i=0
    EhopDNN = []
    while i < len(preDNN):
        TESTE = preDNN[i:len(_R)+i]
        EhopDNN.append(TESTE)
        i+=len(_R)
    
    legendaDNN = ["1 hop(s) NN", "2 hop(s) NN", "4 hop(s) NN", "8 hop(s) NN"]
    legendaMCE = ["1 hop(s) ECM", "2 hop(s) ECM", "4 hop(s) ECM", "8 hop(s) ECM"]
	
    plotGraphDNNvs(legendaMCE, EhopMCE, legendaDNN, EhopDNN, _R, typeG)

#
def plotGraphDNNvs(legendaMCE, EhopMCE, legendaDNN, EhopDNN, _R, typeG):
        
    x = _R
    fonte = 20
    linewidth = 3
    loclegend = 'upper right'
    
    #LABEL
    plt.ylabel('Eihop(dBmJ/bit)', fontsize=fonte)
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
    
    if(typeG != 'R(bps)'):
        loclegend = 'lower left'
        
    #LEGENDAS        
    legenda = []
    legenda.extend(legendaMCE)
    legenda.extend(legendaDNN)
    plt.legend(legenda, ncol=2, loc=loclegend, prop={'size': 12})#upper - lower
    
    #PLOT
    
    plt.savefig('./GRAFICOFINAL.png', format='png', dpi=1000)
    plt.show()
    