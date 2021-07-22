# -*- coding: utf-8 -*-
#print('Omega: \u03A9')     print('Delta: \u0394')      print('sigma: \u03C3')      
#print('mu: \u03BC')        print('epsilon: \u03B5')    print('degree: \u00B0')

import math as mt
from os import path
import pandas as pd

from classFaltasTrifasicas import faltasTrifasicas

sistema = faltasTrifasicas() # Objeto da classe



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~ CAPTURA DADOS DO .CSV EXPORTADO PELO PW ~~~~~~~~~~
barras_df = pd.read_csv('Caminho do seus dados de barras .CSV') # Note que também é possível incluir os dados diretamente através do método .setBarras(numero, tensao, angulo)
barras_df = pd.DataFrame(barras_df) # IMPORTANTE: EXCLUIR O CAMPO 'Bus' DO .CSV! (célula A1) ( O .csv deve ser column header, separado por virgula)

# Insere as barras do sistema
for i in range(len(barras_df)):
    num = int(barras_df['Number'][i])         # Número da barra
    ten = float(barras_df['PU Volt'][i])      # Tensão na barra
    ang = float(barras_df['Angle (Deg)'][i])  # Ângulo da Tensão
    
    sistema.setBarras(num, ten, ang)


# Caputra dados das conexões
conexoes_df = pd.read_csv('Caminho do seus dados de conexões entre as barras .csv') # Note que também é possível incluir os dados diretamente através do método .setConexao(barra de, barra para, impedancia)
conexoes_df = pd.DataFrame(conexoes_df) # IMPORTANTE: EXCLUIR O CAMPO 'Branch' DO .CSV! (célula A1)

for i in range(len(conexoes_df)):
    barra_de = int(conexoes_df['From Number'][i])                  # Barra de
    barra_para = int(conexoes_df['To Number'][i])                  # Barra para
    z_linha = complex(conexoes_df['R'][i], conexoes_df['X'][i])  # Impedância da linha de transmissão

    sistema.setConexao(barra_de, barra_para, impedancia = z_linha)

# Shunt nas barras (não aparece no conexoesSistema.csv)
sistema.shunt = [0.045j, 0j, 0.0225j, 0j, 0j] 

sistema.printBarras()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~ Zbus DO SISTEMA  ~~~~~~~~~~~~~~
sistema.Zbus(printDados = True) # calculo da Zbarra

sistema.bus2Gnd(toBus = 2, Zb=complex(0.2067, 0.0775), printDados=True)   # Zb indutivo entre barra existente e terra


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~ CÁLCULO DA FALTA SIMÉTRICA ~~~~~~~~~~
#sistema.faltaTrifasica(bus = 2, printDados = True) # calcula falta trifásica na barra
for i in range(1, len(barras_df)+1):
    sistema.faltaTrifasica(bus = i, printDados = True) # calcula falta trifásica na barra


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~ CONSTRUÇÃO DIRETA DA Zbus ~~~~~~~~~~

#sistema.bus2Gnd(bus = 2, Zb=complex(0.2067, 0.0775), printDados=True)   # Zb entre barra existente e terra

#sistema.newBus2Gnd(Zb = complex(0, 1.2), printDados=True)               # Zb entre nova barra e referência

#sistema.newBus2Bus(toBus = 2, Zb = complex(0, 2.3)printDados=True)      # Zb entre nova barra em barra existente

#sistema.bus2Gnd(toBus=1, Zb=complex(0, 5), printDados=True)             # Zb entre barra existente e terra

#sistema.bus2Bus(bus=1, toBus=5, Zb=complex(0, 1), printDados=True)


