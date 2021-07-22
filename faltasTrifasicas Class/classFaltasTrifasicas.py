################################################################################################
# @title: Faltas Trifásicas                                                                    #
# @description: Código para realizar o cálculo de faltas simétricas em determinada barra. Além #
#               disso, manipula a matriz Zbus de forma direta.                                 #
# @author: Dihordy Ross - dihordyross@gmailcom                                                 #
# @date:   20-Jul.-2021                                                                        #
# @version: Python 3.7                                                                         #
################################################################################################


import decimal
import cmath as cmt
import math as mt
import numpy as np
from numpy import linalg
from numpy.core.fromnumeric import size
import pandas as pd 


class faltasTrifasicas:
    def __init__(self):
        self.__Sbase = 100e6

        self.__dados     = dict()
        self.__conexao   = dict()
        self.shunt    = list()

        self.__Zbus = list()
        self.__ybus = list()

        self.__V = dict()

        self.__tensaoPlot = dict()
        self.__angPlot = dict()
        
        

    # ===================================================================== #
    def setBarras(self, barra, tensao, ang):
        """
        Método que define as barras do sistema

        :param barra: número da barra
        :param tensao: nível de tensão na barra dado em p.u.
        :param ang: ângulo da tensão na barra dado em graus
        """

        self.__dados[barra]      = {'tensao': tensao,'ang': mt.radians(ang)}
        self.__tensaoPlot[barra] = [tensao]
        self.__angPlot[barra]    = [ang]
                                                

    # ===================================================================== #
    def setConexao(self, barra1, barra2, impedancia=None, admitancia=None):
        """
        Método que define as conexões entre as barras do sistema

        :param barra: número da barra
        :param tensao: nível de tensão na barra
        :param ang: ângulo da tensão na barra
        """

        if impedancia is None:
            impedancia = 1 / admitancia
        elif admitancia is None:
            admitancia = 1 / impedancia
        else:
            return "Erro!"

        self.__conexao[(barra1, barra2)] = {'Impedancia': impedancia, 'Admitancia': admitancia, }  # barra__de -> barra__para
        

    # ===================================================================== #
    def Ybus(self, printDados=False):
        """
        Método que calcula a matriz Ybus do sistema

        :param printDados: True para apresentar valores
        """

        self.__ybus = np.ones((len(self.__dados), len(self.__dados)), dtype = complex)

        for i in range(len(self.__ybus)):
            lin = []

            for j in range(len(self.__ybus)):
                if i == j: # diagonal principal preenchida com zeros
                    lin.append(0)
                else: # diagonal secundaria
                    if self.__conexao.__contains__(tuple([i + 1, j + 1])):
                        lin.append(-self.__conexao.get(tuple([i + 1, j + 1]))['Admitancia'])

                    elif self.__conexao.__contains__(tuple([j + 1, i + 1])):
                        lin.append(-self.__conexao.get(tuple([j + 1, i + 1]))['Admitancia'])
                    else:
                        lin.append(0)

            for j in range(len(self.__ybus)):
                if i==j: 
                    lin[j] = -1 * sum(lin)

            self.__ybus[i] = lin # Ybus

        # inclui elementos shunt
        for i in range(0, len(self.__ybus)):
            try:
                #print(1/self.shunt[i])
                self.__ybus[i][i] = self.__ybus[i][i] + 1/self.shunt[i]
            except: pass
        
        self.__ybus = np.round(self.__ybus, decimals=4)

        if printDados:
            print('\n Ybus: ')
            print(self.__ybus)
        print('\n___________________________________________________________________')  
        


    # ===================================================================== #
    def Zbus(self, printDados=False):
        """
        Método que calcula a matriz Zbus inicial do sistema

        :param printDados: True para apresentar valores
        """

        self.Ybus(printDados=False)
        self.__Zbus = np.linalg.pinv(self.__ybus)

        if printDados:
            self.__printZbus()
        print('\n___________________________________________________________________')  


    # ===================================================================== #
    def faltaTrifasica(self, bus=int, printDados=False):
        """
        Método que calcula uma falta simétrica em uma barra do sistema

        :param bus: barra onde ocorre a falta
        :param printDados: True para apresentar valores
        """

        self.tensaoSistema(printDados=False)

        print("\n\t ~~~~~~ FALTA SIMÉTRICA NA BARRA ", bus , "~~~~~~")
             
        bus = bus-1 # corrige o indice para matriz

        # corrente de falta I''
        Ifalta = self.__V.get(bus+1)/self.__Zbus[bus][bus]
            
        # tensões pós falta
        Vfalta = dict()
        for i in range(1, len(self.__V)+1):
            Vfalta[i] = self.__V.get(i) -  (self.__Zbus[i-1][bus] * Ifalta)


        # corrente de falta nos ramos
        Iramos = dict()
        self.__ybus = np.linalg.pinv(self.__Zbus)
        for i in self.__dados:
            for j in self.__dados:
                if i == j:
                    continue
                else: 
                    Iramos[(i, j)] = (Vfalta.get(i) - Vfalta.get(j)) * self.__ybus[i - 1][j - 1] # Ikm = (Vk-Vm)*Ykm
        
        
        # Apresenta resultados
        if printDados:
            print("\n Tensão nas barras do sistema pós falta: \n")
            for i in range(1, len(Vfalta)+1):
                print("Barra ", i, "\t ", np.round(abs(Vfalta[i]), decimals=4), "/_", np.round((np.degrees(np.angle(Vfalta[i]))), decimals=4), "° ")

            print("\n Corrente de falta trifásica:  \n")
            print(np.round(abs(Ifalta), decimals=4), " /_ ", np.round((np.degrees(np.angle(Ifalta))), decimals=4), "° \n")

            # print correntes diferentes de zero (onde há conexão entre barras)
            for i in self.__dados:
                for j in self.__dados:
                    if i == j:
                        continue
                    else:
                        if Iramos[(i, j)] != 0+0j:
                            print("\n Corrente de falta no ramo: ", i ,",", j ,"\n") 
                            print(np.round(abs(Iramos[(i, j)]), decimals=4), " /_ ", np.round((180/mt.pi)*(np.angle(Iramos[(i, j)])), decimals=4), "° \n")
            print("___________________________________________________________________")

     

    # ===================================================================== #
    def tensaoSistema(self, printDados=False):
        """
        Método que define as tensões nas barras do sistema

        :param printDados: True para apresentar valores
        """
        self.__V = dict()
        for i in self.__dados:
            self.__V[i] = cmt.rect( self.__dados.get(i)['tensao'],
                                    self.__dados.get(i)['ang'])
        if printDados:
            self.__printTensao()


    # ===================================================================== #
    def newBus2Gnd(self, Zb=complex, printDados=False):        
        """
        Método que inclui na Zbus uma impedância entre uma nova barra e a barra de referência

        :param Zb: impedância a ser conectada
        :param printDados: True para apresentar valores
        """

        print("\n Inserindo barra nova ", len(self.__Zbus)+1, " para referência \n")
        aux = np.zeros(shape=(len(self.__Zbus)+1, len(self.__Zbus)+1), dtype=complex)

        for i in range(len(self.__Zbus)):
            for j in range(len(self.__Zbus)):
                aux[i][j] = self.__Zbus[i][j]

        # insere ultimo elemento
        aux[len(aux)-1, len(aux)-1] = Zb
        self.__Zbus = aux
        
        if printDados:
            self.__printZbus()
        print('\n___________________________________________________________________\n')  

        
    # ===================================================================== #
    def newBus2Bus(self, toBus=int, Zb=complex, printDados=False):     
        """
        Método que inclui na Zbus uma impedância entre uma nova barra e uma barra existente

        :param toBus: barra para conexão da impedância (existente)
        :param Zb: impedância a ser conectada
        :param printDados: True para apresentar valores
        """

        print("\n Inserindo barra nova ", len(self.__Zbus)+1, " para barra existente ", toBus ,"\n")
        aux = np.zeros(shape=(len(self.__Zbus)+1, len(self.__Zbus)+1), dtype=complex)
        toBus = toBus - 1  #corrige indice da barra na matriz

        for i in range(len(self.__Zbus)):
            for j in range(len(self.__Zbus)):
                aux[i][j] = self.__Zbus[i][j]
               
        # copia coluna toBus
        for j in range(len(self.__Zbus)):
            aux[len(aux)-1][j] = self.__Zbus[toBus][j]

        # copia linha toBus
        for i in range(len(self.__Zbus)):
            aux[i][len(aux)-1] = self.__Zbus[i][toBus]

        # insere ultimo elemento
        aux[len(aux)-1, len(aux)-1] = self.__Zbus[toBus][toBus] + Zb
        self.__Zbus = aux

        
        if printDados:
            self.__printZbus()
        print('\n___________________________________________________________________\n')  


    # ===================================================================== #
    def bus2Gnd(self, toBus=int, Zb=complex, printDados=False):
        """
        Método que insere na Zbus uma impedância entre uma barra existente e a barra de referência

        :param toBus: barra para conexão da impedância (existente)
        :param Zb: impedância a ser conectada
        :param printDados: True para apresentar valores
        """

        print("\n Inserindo entre barra existente ", toBus, " para referência \n")

        aux = np.zeros(shape=(len(self.__Zbus)+1, len(self.__Zbus)+1), dtype=complex)
        toBus = toBus - 1  #corrige indice da barra na matriz

        for i in range(len(self.__Zbus)):
            for j in range(len(self.__Zbus)):
                aux[i][j] = self.__Zbus[i][j]
               
        # copia coluna toBus
        for j in range(len(self.__Zbus)):
            aux[len(aux)-1][j] = self.__Zbus[toBus][j]

        # copia linha toBus
        for i in range(len(self.__Zbus)):
            aux[i][len(aux)-1] = self.__Zbus[i][toBus]

        # insere ultimo elemento
        aux[len(aux)-1, len(aux)-1] = self.__Zbus[toBus][toBus] + Zb
        self.__Zbus = aux

        if printDados:
            self.__printZbus()
        print('\n___________________________________________________________________\n')  


        # Redução de Kron
        self.kron(k=len(self.__Zbus), printDados=True)



    # ===================================================================== #
    def bus2Bus(self, bus=int, toBus=int, Zb=complex, printDados=False):
        """
        Método que insere na Zbus uma impedância entre duas barras existentes

        :param bus: barra existente 'j' para conexão
        :param toBus: barra existente 'k' para conexão
        :param Zb: impedância a ser conectada
        :param printDados: True para apresentar valores
        """

        print("\n Inserindo entre barra existente ", bus, " e barra existente ", toBus ,"\n")
        bus = bus - 1  #corrigindo indices para matriz
        toBus = toBus - 1

        aux = np.zeros(shape=(len(self.__Zbus)+1, len(self.__Zbus)+1), dtype=complex)

        for i in range(len(self.__Zbus)):
            for j in range(len(self.__Zbus)):
                aux[i][j] = self.__Zbus[i][j]
               
        # Coluna j - Coluna k
        for j in range(len(self.__Zbus)):
            aux[j][len(aux)-1] = self.__Zbus[j][bus] - self.__Zbus[j][toBus]

        # Linha j - Linha k
        for i in range(len(self.__Zbus)):
            aux[len(aux)-1][i] = self.__Zbus[bus][i] - self.__Zbus[toBus][i]

        # insere ultimo elemento
        aux[len(aux)-1, len(aux)-1] = self.__Zbus[bus][bus] + self.__Zbus[toBus][toBus] - 2*self.__Zbus[bus][toBus] + Zb

        if printDados:
            self.__printZbus()
        print('\n___________________________________________________________________\n') 

        self.kron(k=len(self.__Zbus), printDados=True)


    # ===================================================================== #
    def kron(self, k=int, printDados=False):      
        """
        Método que realiza a redução de Kron, removendo barras do sistema

        :param k: barra existente 'k' para ser removida
        :param Zb: impedância a ser conectada
        """
        

        k = k-1 #corrige o indice para a matriz
        Znovo = np.zeros(shape=(len(self.__Zbus)-1, len(self.__Zbus)-1), dtype=complex)
        
        for i in range(0, len(self.__Zbus)-1):
            for j in range(0, len(self.__Zbus)-1):
                Znovo[i][j] = self.__Zbus[i][j] - ((self.__Zbus[i][k]*self.__Zbus[k][j]) / self.__Zbus[k][k])

        self.__Zbus = Znovo
        
        if printDados:
            print('\n Redução de Kron para barra ', k, '\n')
            #print("Ybus pós Kron \n", np.round(np.linalg.pinv(self.__Zbus), decimals=4))
            self.__printZbus()
        print('\n___________________________________________________________________\n')  



    # ===================================================================== #
    # Métodos que apresentam valores
    def __printZbus(self):
        print("Zbus: \n", np.round((self.__Zbus), decimals=4))
        print(size(self.__Zbus, axis=0), "x", size(self.__Zbus, axis=1))


    def __printTensao(self):
        print('___________________________________________________________________')  
        print("Tensão nas barras do sistema: \n")
        for i in self.__V:
            print('BARRA \t', i, '\t  TENSÃO = \t', round(abs(self.__V.get(i)), 4), ' /_ ', round(np.degrees(np.angle(self.__V.get(i))), 4), '° \t [pu]')

        print('___________________________________________________________________')





    def printBarras(self):
                    
        print('___________________________________________________________________')  
        print('\n Barras: ')

        for i in self.__dados: print(self.__dados[i])
        print('___________________________________________________________________')  
        
        print('\n Conexões: ')
        for i in self.__conexao: print('Conexão ', i, '\t', self.__conexao[i])
        print('___________________________________________________________________')  


