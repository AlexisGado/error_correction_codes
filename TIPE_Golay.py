import pylab as pl
import time
import random


tot = time.time()


def modif(m,i):
    m[i] = (m[i]+1)%2

# CODE DE GOLAY 

A = pl.matrix([[1,1,1,1,1,1,1,1,1,1,1],
        [1,1,0,1,1,1,0,0,0,1,0],
        [1,0,1,1,1,0,0,0,1,0,1],
        [0,1,1,1,0,0,0,1,0,1,1],
        [1,1,1,0,0,0,1,0,1,1,0],
        [1,1,0,0,0,1,0,1,1,0,1],
        [1,0,0,0,1,0,1,1,0,1,1],
        [0,0,0,1,0,1,1,0,1,1,1],
        [0,0,1,0,1,1,0,1,1,1,0],
        [0,1,0,1,1,0,1,1,1,0,0],
        [1,0,1,1,0,1,1,1,0,0,0],
        [0,1,1,0,1,1,1,0,0,0,1]])
        

gene = pl.block([pl.eye(12),A])
contro = pl.block([[A],[pl.eye(11)]])

print(gene)
# print("")
# print(contro)
# print("")

def code(m):
    res = pl.dot(m,gene)
    for i in range(len(res)):
        res[i]= res[i]%2
    return pl.ravel(res)
    


def syndrome(m):
    res = pl.dot(m,contro)
    for i in range(len(res)):
        res[i]= res[i]%2
    return pl.ravel(res)
    

    

def creaTable() :
    LUTs = []
    LUTs.append([pl.zeros((23)),pl.zeros((11))])
    for i in range(23):
        l=pl.zeros((23))
        modif(l,i)
        LUTs.append([l,syndrome(l)])
        
    for i in range(23):
        for j in range(i+1,23):
            l=pl.zeros((23))
            modif(l,i)
            modif(l,j)
            LUTs.append([l,syndrome(l)])
        
    for i in range(23):
        for j in range(i+1,23):
            for k in range(j+1,23):
                l=pl.zeros((23))
                modif(l,i)
                modif(l,j)
                modif(l,k)
                LUTs.append([l,syndrome(l)])
    return pl.array(LUTs)
    
    
t = time.time()
Table = creaTable()
print("temps table golay :")
print(time.time()-t)
print("")

def decode(m):
    synd = syndrome(m)
    for e in Table:
        
        if (synd == e[1]).all():
            err = e[0]
            corr  = err + m 
            for i in range(len(corr)):
                corr[i]= corr[i]%2
            return corr[:12]
    print("erreur decodage golay")


# CODE DE HAMMING

B = pl.matrix([[0,1,1],[1,0,1],[1,1,0],[1,1,1]])

geneh = pl.block([pl.eye(4),B])
controh = pl.block([[B],[pl.eye(3)]])


def codeh(m):
    res = pl.dot(m,geneh)
    for i in range(len(res)):
        res[i]= res[i]%2
    return pl.ravel(res)
    


def syndromeh(m):
    res = pl.dot(m,controh)
    for i in range(len(res)):
        res[i]= res[i]%2
    return pl.ravel(res)
    
    

def creaTableh() :
    LUTs = []
    LUTs.append([pl.zeros((7)),pl.zeros((3))])
    for i in range(7):
        l=pl.zeros((7))
        modif(l,i)
        LUTs.append([l,syndromeh(l)])
        
    return pl.array(LUTs)
    
    
t = time.time()
Tableh = creaTableh()

print("temps table hamming :")
print(time.time()-t)
print("")

def decodeh(m):
    synd = syndromeh(m)
    for e in Tableh:
        
        if (synd == e[1]).all():
            err = e[0]
            corr  = err + m 
            for i in range(len(corr)):
                corr[i]= corr[i]%2
            return corr[:4]
    print("erreur decodage hamming")







# IMAGES / MODIFICATIONS 



def binim(im):
    p,q = im.shape[0],im.shape[1]
    nveau = pl.zeros((p,q))
    for i in range(p):
        for j in range(q):
            
            if pl.sum(im[i,j])> 1.5:
                nveau[i,j] = 1
            else :
                nveau[i,j] = 0
    return nveau
    
def imaff(bin):
    p,q = bin.shape[0],bin.shape[1]
    nveau = pl.zeros((p,q,3))
    for i in range(p):
        for j in range(q):
            a = bin[i,j]
            nveau[i,j]= pl.array([a,a,a])
            
    return nveau
    
    
pingu = pl.imread('/Users/alexgado/Documents/INFO/IPT/pinguin.jpg')
auguste = pl.imread('/Users/alexgado/Documents/INFO/IPT/auguste.JPG')    
axiome = pl.imread('/Users/alexgado/Documents/INFO/IPT/axiome.png')
ellipse = pl.imread('/Users/alexgado/Documents/INFO/IPT/ellipse.png')
panda = pl.imread('/Users/alexgado/Documents/INFO/IPT/panda.png')



def modifvect(m,p):
    for i in range(len(m)):
        r = random.random()
        if r < p:
            m[i] = (m[i]+1)%2 
            

def modifim(im,pr):
    bin = binim(im)
    (p,q) = bin.shape
    ligne = pl.ravel(pl.reshape(bin, (1,p*q)))
    modifvect(ligne,pr)
    finalbin = pl.reshape(ligne,(p,q))
    return imaff(finalbin)



# TESTS GOLAY



def modifimcorr(im,pr):
    bin = binim(im)
    (p,q) = bin.shape
    ligne = pl.ravel(pl.reshape(bin, (1,p*q)))
    nveau = []
    for i in range(0,int(p*q)-11,12):
        
        bloc = ligne[i:i+12]
        codeb = code(bloc)

        modifvect(codeb,pr)       
        deco = decode(codeb)
        nveau.append(deco)
    nveau.append(ligne[p*q-1:p*q-1-p*q%12:-1])
    
    lis = []
    for m in nveau:
        for i in range(len(m)):
            lis.append(m[i])
    finalbin = pl.reshape(pl.array(lis),(p,q))
    return imaff(finalbin)
    


    
# TESTS HAMMING

def modifimcorrh(im,pr):
    bin = binim(im)
    (p,q) = bin.shape
    ligne = pl.ravel(pl.reshape(bin, (1,p*q)))
    nveau = []
    for i in range(0,int(p*q)-3,4):
        
        bloc = ligne[i:i+4]
        codeb = codeh(bloc)

        modifvect(codeb,pr)       
        deco = decodeh(codeb)
        nveau.append(deco)
    nveau.append(ligne[p*q-1:p*q-1-p*q%4:-1])
    
    lis = []
    for m in nveau:
        for i in range(len(m)):
            lis.append(m[i])
    finalbin = pl.reshape(pl.array(lis),(p,q))
    return imaff(finalbin)



# COMPARAISON
    
def compare(im,pr):
    pl.figure()
    pl.imshow(imaff(binim(im)))
    pl.show()
    pl.figure()
    pl.imshow(modifim(im,pr))
    pl.show()
    pl.figure()
    pl.imshow(modifimcorr(im,pr))
    pl.show()
    pl.figure()
    pl.imshow(modifimcorrh(im,pr))
    pl.show()



compare(ellipse,0.1)

print("temps total :")
print(time.time()-tot)


    



























