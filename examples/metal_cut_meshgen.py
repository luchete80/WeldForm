#https://www.tutorialsteacher.com/python/create-ui-using-tkinter-in-python

from tkinter import *
from tkinter.ttk import Combobox
import numpy as np
import math

window=Tk()
linea_g=10

filelbl = Label(window, text="Archivo", width=15,justify=LEFT)
filelbl.grid(column=1, row=0)	
textField = Entry(window, width=15)
textField.grid(column=2, row=0)
textField.insert(0,"test.txt")
#textField.pack(fill=NONE, side=TOP)

lbl = Label(window, text="Largo [mm]")
lbl.grid(column=1, row=1)
largoField = Entry(window, width=10)
largoField.insert(0, "100.")
largoField.grid(column=2, row=1)

lbl = Label(window, text="Ancho [mm]")
lbl.grid(column=1, row=2)
anchoField = Entry(window, width=10)
anchoField.insert(0, "12.")
anchoField.grid(column=2, row=2)
#largoField.pack(fill=NONE, side=TOP)

lbl = Label(window, text="Paso [mm]")
lbl.grid(column=1, row=3)
pasoField = Entry(window, width=10)
pasoField.insert(0, "4.0")
pasoField.grid(column=2, row=3)

lbl = Label(window, text="Pasadas [mm]")
lbl.grid(column=1, row=4)
npasField = Entry(window, width=10)
npasField.insert(0, "10")
npasField.grid(column=2, row=4)

lbl = Label(window, text="Paso Z [mm]")
lbl.grid(column=1, row=5)
pasoZField = Entry(window, width=10)
pasoZField.insert(0, "2.0")
pasoZField.grid(column=2, row=5)

lbl = Label(window, text="Long lenta [mm]")
lbl.grid(column=1, row=6)
pasoZField = Entry(window, width=10)
pasoZField.insert(0, "10.0")
pasoZField.grid(column=2, row=6)

lbl = Label(window, text="V Lenta [%]")
lbl.grid(column=1, row=7)
pasoZField = Entry(window, width=10)
pasoZField.insert(0, "50")
pasoZField.grid(column=2, row=7)


#INPUTS 
L = 1
dx = 0.1


node = [ [0]*2 for i in range(10)]


alpha = 15.0 * np.pi/180.0 
nside = L/dx
x = L
y = L*math.cos(alpha)

print("alpaa: ", alpha)
for i in range(int(nside)):
    x -=dx * math.cos(alpha)
    y-= dx * math.sin(alpha)
    node[i][0] = x 
    node[i][1] = y - dx * math.cos(alpha)
print("Node ",node[:][:])
    
    

#PRINT ELEMENT 
#https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/element.pdf


#********** ARCHIVO *************

# added only to define required global variable "txt"
class Txt(object):
    def SetValue(data): pass
    def GetValue(): pass
txt = Txt()
####

def load():
    with open(textField.get()) as file:
        f= open(file.read(),"w+")
		
def lin2str(l):
    ret=""
    if (l < 100):
        ret="N000"+str(l)
    if (l < 100):
        ret="N00"+str(l)
    if (l < 1000):
        ret="N0"+str(l)
    else:
        ret="N"+str(l)
    l=l+5
    return ret
    
    

# Muy bueno para devolver variables
# https://stackoverflow.com/questions/15286401/print-multiple-arguments-in-python
def save(lin):
    lin=10
    largo=float(largoField.get())
    paso =float(pasoField.get())
    pasoZ=float(pasoZField.get())
    ancho=float(anchoField.get())
    npasadas=int(npasField.get())
    npasos=int(largo/paso)
    print("Largo: %.2f, Ancho: %.2f" %(largo,ancho))
    print("Pasos: %d" %npasos)
    #with open(textField.get(), 'w') as file:
	#with open(textField.get(), 'w') as file:
    f= open(textField.get(),"w+")

    z=0.5	
    #generate mesh
    for p in range(npasadas):
	    #Veo para que lado va la pasada   
        if p%2==0: #Pasada par
    	    fs=1.0
        else:
            fs=-1.0
        x=-fs*largo/2.0	
        y=-fs*ancho/2.0
        f.write(lin2str(lin)+" X%.4f Y%.4f\n" % (x,y))
        lin=lin+10	
        f.write(lin2str(lin)+" Z%.4f \n" % (z))
        lin=lin+10
        f.write(lin2str(lin)+" M11P13 \n")
        for i in range(npasos):
            y=ancho/2.
            f.write(lin2str(lin)+" Y%.4f \n" % (y))
            lin=lin+10
            x=x+fs*paso/2.
            f.write(lin2str(lin)+" X%.4f \n" % (x))
            lin=lin+10
            y=-ancho/2.
            f.write(lin2str(lin)+" Y%.4f \n" % (y))
            lin=lin+10
            x=x+fs*paso/2.
            f.write(lin2str(lin)+" X%.4f \n" % (x))	
            lin=lin+10
			# create text field
        f.write(lin2str(lin)+"( ******************* FIN DE PASADA ***********************)\n")
        lin=lin+10
        f.write(lin2str(lin)+" M10P13 \n")
        lin=lin+10
        f.write(lin2str(lin)+" G4 P60000\n")
        lin=lin+10
        z=z+pasoZ
    f.close()

		
#Si no se coloca lambda no funciona
b = Button(window, text="Generar", width=10, command=lambda:save(linea_g))
b.grid(column=3, row=10)
#b.pack()

window.title('Generador Codigo G')
window.geometry("400x200+10+10")
window.mainloop()