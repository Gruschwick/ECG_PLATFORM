import numpy as np
from scipy import signal
import math
import wfdb

def lynfilt (x, fs, ns):
    
    ns = np.size(x)
    
    Bd = [[1,2,0,-2,-1]]
    Bd = (fs/8) * Bd
    Ad = 1
    Der = signal.lfilter(Bd, Ad, x)
    Der[2] = fs * (x[1] - x[0])
    Der[3] = (fs/4) * (2*x[2]*x[0])
    Td = 2
    
    T = Td
    Der[0:ns-T] = Der[T:ns]
    Der[ns-T:ns] = np.zeros(T,1)
    
    rmax = max(abs(Der[0:2*fs]))
    if rmax == 0:
        rmax =1
    Der = np.divide(10*Der,rmax)
    
    Fpa = 1
    mpa = round(fs/Fpa)
    Bpa = np.zeros(1,2*(mpa+1))
    Bpa[0] = 1
    Bpa[mpa] = -2
    Bpa[2*mpa] = 1
    Apa = [[1,-2,1]]
    Xpa = signal.lfilter(Bpa,Apa,x)
    Tpa = (mpa - 1)
    
    T=Tpa
    Xpa[0:ns-T] = (mpa)*(mpa) * x[0:ns-T] - Xpa[T:ns]
    Xpa[ns-T:ns] = np.zeros(T,1)
    
    rmax = max(abs(Xpa[0:2*fs]))
    if rmax == 0:
        rmax =1
    Xpa = np.divide(10*Xpa,rmax)
    
    Fpb = 60
    mpb = round(fs/Fpb)
    Bpb = np.zeros(1,2*mpb+1)
    Bpb[0] = 1
    Bpb[mpb] = -2
    Bpb[2*mpb] = 1
    Apb = [[1,-2,1]]
    Xpb = signal.lfilter(Bpb,Apb,x)
    Tpb = (mpb - 1)
    
    T=Tpb + 1
    Xpb[0:ns-T] = Xpa[T:ns]
    Xpb[ns-T:ns] = np.zeros(T,1)
    
    rmax = max(abs(Xpb[0:2*fs]))
    if rmax == 0:
        rmax =1
    Xpb = np.divide(10*Xpb,rmax)
    
    Bd = [[1,2,0,-2,-1]]
    Bd = (fs/8)*Bd
    Ad = 1
    D = signal.lfilter(Bd,Ad,Xpb)
    D[2] = fs * (Xpb[1]-Xpb[0])
    D[3] = (fs/4) * (2*Xpb[2] - 2*Xpb[1])
    Td = 2
    
    T = Td
    D[0:ns-T] = D[T:ns]
    D[ns-T:ns] = np.zeros(T,1)

    rmax = max(abs(D[0:2*fs]))
    if rmax == 0:
        rmax =1
    D = np.divide(10*D,rmax)
    
    Fpf = 40
    mpf = round(fs/Fpf)
    Bpf = np.zeros(1,2*mpb+1)
    Bpf[0] = 1
    Bpf[mpf] = -2
    Bpf[2*mpf] = 1
    Apf = [[1,-2,1]]
    Xpf = signal.lfilter(Bpf,Apf,x)
    Tpf = (mpf - 1)
    
    T=Tpf + 1
    Xpf[0:ns-T] = Xpf[T:ns]
    Xpf[ns-T:ns] = np.zeros(T,1)
    
    rmax = max(abs(Xpf[0:2*fs]))
    if rmax == 0:
        rmax =1
    Xpf = np.divide(10*Xpf,rmax)
    
    F = signal.lfilter(Bd,Ad,Xpf)
    F[2] = fs * (Xpf[1]-Xpf[0])
    F[3] = (fs/4) * (2*Xpf[2] - 2*Xpf[1])
    Td = 2
    
    T = Td
    F[0:ns-T] = F[T:ns]
    F[ns-T:ns] = np.zeros(T,1)

    rmax = max(abs(F[0:2*fs]))
    if rmax == 0:
        rmax =1
    F = np.divide(10*F,rmax)  
    
    return Xpa,Xpb,D,F,Der

def lymfilt2(X,fs,ns):
    
    ns = np.size(X)
    
    Bd = [[1,2,0,-2,-1]]
    Bd = (fs/8) * Bd
    Ad = 1
    Der = signal.lfilter(Bd, Ad, X)
    Der[2] = fs * (X[1] - X[0])
    Der[3] = (fs/4) * (2*X[2] - 2*[0])
    Td = 2
    
    T = Td
    Der[0:ns-T] = Der[T:ns]
    Der[ns-T:ns] = np.zeros(T,1)
    
    if np.size(Der) >= 2*fs:
        
        rmax = max(abs(Der[0:2*fs]))
    
    else:
        rmax = max(abs(Der))

    if rmax == 0:
        rmax =1
    Der = np.divide(10*Der,rmax)
    
    Fpa = 1
    mpa = round(fs/Fpa)
    Bpa = np.zeros(1,2*mpa+1)
    Bpa[0] = 1
    Bpa[mpa] = -2
    Bpa[2*mpa] = 1
    Apa = [[1,-2,1]]
    Xpa = signal.lfilter(Bpa,Apa,X)
    Tpa = (mpa - 1)
    
    T=Tpa
    Xpa[0:ns-T] = (mpa)*(mpa) * X[0:ns-T] - Xpa[T:ns]
    Xpa[ns-T:ns] = np.zeros(T,1)

    
    if np.size(Xpa) >= 2*fs:
        
        rmax = max(abs(Xpa[0:2*fs]))
    
    else:
        rmax = max(abs(Xpa))

    if rmax == 0:
        rmax =1
    Xpa = np.divide(10*Xpa,rmax)
    
    Fpb = 45
    mpb = round(fs/Fpb)
    Bpb = np.zeros(1,2*mpb+1)
    Bpb[0] = 1
    Bpb[mpb] = -2
    Bpb[2*mpb] = 1
    Apb = [[1,-2,1]]
    Xpb = signal.lfilter(Bpb,Apb,X)
    Tpb = (mpb - 1)
    
    T=Tpb + 1
    Xpb[0:ns-T] = Xpa[T:ns]
    Xpb[ns-T:ns] = np.zeros(T,1)
    
    if np.size(Xpb) >= 2*fs:
        
        rmax = max(abs(Xpb[0:2*fs]))
    
    else:
        rmax = max(abs(Xpb))
    if rmax == 0:
        rmax =1
    Xpb = np.divide(10*Xpb,rmax)
    
    Bd = [[1,2,0,-2,-1]]
    Bd = (fs/8)*Bd
    Ad = 1
    D = signal.lfilter(Bd,Ad,Xpb)
    D[2] = fs * (Xpb[1]-Xpb[0])
    D[3] = (fs/4) * (2*Xpb[2] - 2*Xpb[1])
    Td = 2
    
    T = Td
    D[0:ns-T] = D[T:ns]
    D[ns-T:ns] = np.zeros(T,1)

    if np.size(D) >= 2*fs:
        
        rmax = max(abs(D[0:2*fs]))
    
    else:
        rmax = max(abs(D))
    if rmax == 0:
        rmax =1
    D = np.divide(10*D,rmax)
    
    Fpf = 40
    mpf = round(fs/Fpf)
    Bpf = np.zeros(1,2*mpb+1)
    Bpf[0] = 1
    Bpf[mpf] = -2
    Bpf[2*mpf] = 1
    Apf = [[1,-2,1]]
    Xpf = signal.lfilter(Bpf,Apf,X)
    Tpf = (mpf - 1)
    
    T=Tpf + 1
    Xpf[0:ns-T] = Xpf[T:ns]
    Xpf[ns-T:ns] = np.zeros(T,1)
    
    if np.size(Xpf) >= 2*fs:
        
        rmax = max(abs(Xpf[0:2*fs]))
    
    else:
        rmax = max(abs(Xpf))
    if rmax == 0:
        rmax =1
    Xpf = np.divide(10*Xpf,rmax)
    
    F = signal.lfilter(Bd,Ad,Xpf)
    F[2] = fs * (Xpf[1]-Xpf[0])
    F[3] = (fs/4) * (2*Xpf[2] - 2*Xpf[1])
    Td = 2
    
    T = Td
    F[0:ns-T] = F[T:ns]
    F[ns-T:ns] = np.zeros(T,1)

    if np.size(F) >= 2*fs:
        
        rmax = max(abs(F[0:2*fs]))
    
    else:
        rmax = max(abs(F))
    if rmax == 0:
        rmax =1
    F = np.divide(10*F,rmax)  
    
    return Xpa,Xpb,D,F,Der
    
def calcrr (RRm,n,fs,t):
    
    perc= 0.5
    
     
    if n <= 1:
         
        if RRm == 0:
             
            RRm = t[1] - t[0]
             
            if (RRm > 1000e-3* fs) or (RRm < 500e-3 * fs):
                 
                RRm = round(650e-3 * fs)
                 
    else:
         
         if n < 7:
             
             if (t[n]-t[n-1]) > RRm*(1-perc) & (t[n]-t[n-1])<RRm*(1+perc):
                 
                 RRm=(RRm*(n-2)/(n-1))+((t[n]-t[n-1])/(n-1));
         else:
                 
            if (t[n]-t[n-1])>RRm*(1-perc) & (t[n]-t[n-1])<RRm*(1+perc):
                     
                RRm=(RRm*4/5)+((t[n]-t[n-1])/5);
                     
    return RRm

def buscacero(X):
    
    I=np.nonzero(np.sign(X) != np.sign(X(0)))
    
    if not not I:
        
        ncero=I[0]
        ncero = ncero-(abs(X[ncero-1]) < abs(X[ncero]))
        
    else: 
        ncero=[]

    return ncero

def crearumbral(X,umbral):
    

    if not not X:
        
        if X[0] > umbral:
        
            I=np.nonzero(X<umbral)
            if not not I:
                
                iumb=I[0]
                
                if abs(X[iumb-2]-umbral)<abs(X[iumb-1]-umbral):
                    iumb=iumb-1;
                    
            else:
                iumb=[]
            
        else:
            
            I=np.nonzero(X>umbral)
            
            if not not I:
                
                iumb=I[0]
                if abs(X[iumb-2]-umbral)<abs(X[iumb-1]-umbral):
                    iumb=iumb-1;
      
    else:
       iumb=[];
   
    return iumb
    
def zerocross(X):
     
     I = np.nonzero (np.sign(X) != np.sign(X(0)))
     
     if not not I:
         
         ncero  = I[0]
         ncero = ncero = (abs(X[ncero-2]) < abs(X[ncero-1]))
         
     else:
         ncero = []
            
     return ncero
 
def thresholdcross(X, umbral):
    
    if not not X:
        
        if X[0] > umbral:
            
            I = np.nonzero(X < umbral)
            
            if not not I:
                
                iumb = I[0]
            
                if  abs(X[iumb-2] - umbral) < abs (X[iumb-1] - umbral):
                
                    iumb = iumb -1
                
            else:
            
                iumn = []
                
        else:
            
            I = np.nonzero(X > umbral)
            
            if not not I:
                
                iumb = I[0]
                
                if  abs(X[iumb-2] - umbral) < abs (X[iumb-1] - umbral):
                
                    iumb = iumb -1
                
            else:
            
                iumb = []
                
    
    return iumb

def buscaruido(X):
    
    inicio = 1
    iew = np.size(X)
    ruido = 0 
    i = 0
    ymin2 = []
    # imin2 = []
    ymax2 = []
    # imax2 = []
    
    while inicio < iew:
        
        ifinal = min(inicio+5,iew)
        # np.concatenate(ymin2,imin2) = min(X[inicio-1:ifinal])??
        ymin2 = min(X[inicio-1:ifinal])
        ymax2 = max(X[inicio-1:ifinal])
        rudio = (ymax2-ymin2) + ruido
        inicio = ifinal
        i +=1
        
    if i>0:
        ruido = abs(ruido/i)
        
    return ruido

def buscapic1(x,thresh,number,sortem, *args):
    
    varargin = args
    nargin = 4 + len(varargin)
    
    if nargin == 1:
        
        thresh = -math.inf
        number = -1
        sortem = 0
    elif nargin == 2:
        number = -1
        sortem = 0 
    elif nargin == 3:
        sortem = 0
        
#     if (strcmp(sortem,'sort'))
#     sortem=1;
# end

    M = np.size(x)
    #X = x.flatte(1)
    
    if any(np.imag(x.flatten(1))) != 0 :
        x = abs(x)
        
    mask = np.diff(np.sign(np.diff(x.flatten(1))))
    mask = np.hstack([0,mask,0])
    jkl = np.where(abs(mask) > 0 and abs(x) >=thresh)[0]
    if number > 0 and np.size(jkl)>number:
        ii = np.sort(abs(x[jkl]))
        ii = np.flipud(ii)
        jkl = jkl[ii[0:number]]
        jkl = np.sort(jkl)
        
    L = np.size(jkl)
    peaks = []
    peaks[0:L]= x[jkl]
    locs = []
    locs[0:L]= jkl
    
    return peaks, locs

def buscapic2(x,thresh,number,sortem, *args):
    
    varargin = args
    nargin = 4 + len(varargin)
    
    if nargin == 1:
        
        thresh = -math.inf
        number = -1
        sortem = 0
    elif nargin == 2:
        number = -1
        sortem = 0 
    elif nargin == 3:
        sortem = 0
        
#     if (strcmp(sortem,'sort'))
#     sortem=1;
# end

    M = np.size(x)
    #X = x.flatte(1)
    
    if any(np.imag(x.flatten(1))) != 0 :
        x = abs(x)
        
    mask = np.diff(np.sign(np.diff(x.flatten(1))))
    mask = np.hstack([0,mask,0])
    jkl = np.where(abs(mask) > 0 and abs(x) >=thresh)[0]
    njkl = np.zeros(1,number)
    if number > 0 and np.size(jkl)>number:
        ii = np.sort(abs(x[jkl]))
        ii = np.flipud(ii)
        #jkl = jkl[ii[0:number]]
        #jkl = np.sort(jkl)
        njkl[0] = jkl[ii[0]]
        i=2
        l=2
        #temp1 = ii[l]
        #temp2 = jkl[temp1]
        #temp3 = x[temp2]
        #temp1x = ii[l-1]
        #temp2x = jkl[temp1x]
        #temp3x = x[temp2x]
        #temp1y = ii[l+1]
        #temp2y = jkl[temp1y]
        #temp3y = x[temp2y]
        
        while i in range(l, np.size(jkl)-1):
            if (abs(x[jkl[ii[l]]]) >= 0.5*abs(x[jkl[ii[l-1]]])) and (abs(x[jkl[ii[l]]]) >=0.5*abs(x[jkl[ii[l+1]]])):
            #if (abs(temp3) >= 0.5*abs(temp3x)):
                if i<=number:
                    njkl[i-1] = jkl[ii[l]]
                    i+=i
        njkl = np.sort(njkl)
        
    L = np.size(jkl)
    peaks = []
    peaks[0:L]= x[njkl]
    locs = []
    locs[0:L]= njkl
    
    return peaks, locs

def peaksearch(x,thresh,number,sortem, *args):
    
    varargin = args
    nargin = 4 + len(varargin)
    
    if nargin == 1:
        
        thresh = -math.inf
        number = -1
        sortem = 0
    elif nargin == 2:
        number = -1
        sortem = 0 
    elif nargin == 3:
        sortem = 0
        
    M = np.size(x)
    if any(np.imag(x.flatten(1))) != 0 :
        x = abs(x)
        
    mask = np.diff(np.sign(np.diff(x.flatten(1))))
    mask = np.hstack([0,mask,0])
    jkl = np.where(abs(mask) > 1 and abs(x) >=thresh)[0]
    #njkl = np.zeros(1,number)
    if number > 0 and np.size(jkl)>number:
        ii = np.sort(abs(x[jkl]))
        ii = np.flipud(ii)
        jkl = jkl[ii[0:number]]
        jkl = np.sort(jkl)

    L = np.size(jkl)
    peaks = []
    peaks[0:L] = x[jkl]
    locs = []
    locs[0:L] = jkl
    
    return peaks, locs
    
        
def testpic(X,pic,Fs,p):
    
    kpos = round(10e-3*Fs)
    laux = pic
    Xaux = X[pic-kpos-1: pic+kpos]
    if p==1:
        iaux = max(Xaux)
        
    else:
        iaux = min(Xaux)
        
    pic = pic - kpos + iaux -1
    
    return pic

def testpeak(X,pic,Fs,p):
    
    kpos = round(10e-3*Fs)
    laux = pic
    Xaux = X[pic-kpos-1:pic+kpos]
    if p==1:
        iaux = max(Xaux)
    else:
        iaux = min(Xaux)
        
    pic = pic-kpos+iaux-1

def  noiselevel(X):
    
    inicio = 1
    iew = np.size(X)
    mnoise = 0
    i=0
    
    while (inicio < iew):
        ifinal = min(inicio+5,iew)
        ymin2 = min(X[inicio-1:ifinal])
        ymax2 = max(X[inicio-1:ifinal])
        mnoise = (ymax2 - ymin2) + mnoise
        inicio = ifinal
        
        i+=i
        
    if i>0:
        mnoise = abs(mnoise/i)
        
    return mnoise
        
def pbound(n,X,Xpb,F,PKni,Rp,QRS1,prevt,dermax,Fs,Kpb,Kpe):
    
    Pex =1
    P1 = []
    Pp = []
    P2 = []
    
    type = 0
    nofi =1
    bwindp = 200e-3
    ewindp = 30e-3
    iew = QRS1 -round(ewindp*Fs)
    ibw = QRS1 - round(bwindp*Fs)
    
    if ibw<=0:
        ibw = 1
    
    while (nofi==1) and not Pp and (QRS1-ibw)/Fs<300e-3:
        if n==1 or prevt == 0:
            nofi = 0
        elif ibw<prevt:
            ibw=prevt
            nofi=0
        else:
            nofi=1
            
    if ibw>iew:
        Pex = 0
        type =2
    else:
        ymin = min(F[ibw-1:iew])
        imin = min(F[ibw-1:iew])
        imin = imin + ibw -1
        ymax = max(F[ibw-1:iew])
        imax = max(F[ibw-1:iew])
        imax = imax + ibw -1
            
        if (imin <=imax) and (QRS1-imax)/Fs <30e-3:
            iew = QRS1 - round(30e-3*Fs)
            ibw = ibw -round(30e-3*Fs)
            ymin = min(F[ibw-1:iew])
            imin = min(F[ibw-1:iew])
            imin = imin +ibw-1
            ymax = max(F[ibw-1:iew])
            imax = max(F[ibw-1:iew])
            imax = imax +ibw-1
        
    Xaux = Xpb[QRS1-round(15e-3*Fs)-1:QRS1-1]
    base = np.mean(Xaux)
    ecgpbmax = max(abs(Xpb[ibw-1:iew] -base))
        
    if (ecgpbmax<=abs(Xpb[PKni-1]-base)/30) or ((ymax<dermax/(100) and abs(ymin)<dermax/(100)) and (ymax<abs(ymin)/1.5 or ymax>abs(ymin)*1.5)) or  (ymax<0 or ymin>0): 
        Pex=0
        type=2
    elif imin<=imax: 
        type=1
        iaux=imin
        yaux=ymin
        imin=imax
        ymin=ymax
        imax=iaux
        ymax=yaux
      
   #% ---- P wave onset ----
        umbral=(ymax/Kpb)
        Faux=np.flipud(F[0:imax])
        iumb=thresholdcross(Faux,umbral);
        if not iumb:
            iumb=imax
        else:
            iumb=imax-iumb+1
            while Pex==1 and ((QRS1-iumb)/Fs>=240e-3 or iumb<=prevt):
                ibw=ibw+20;
                if ibw>iew-round(20e-3*Fs):
                    Pex=0 
                    type=2;
                else:
                    ymin2 = min(F[ibw-1:iew])
                    imin2 = min(F[ibw-1:iew])
                    imin2=imin2+ibw-1;
                    ymax2 = max(F[ibw-1:iew])
                    imax2 = max(F[ibw-1:iew])
                    imax2=imax2+ibw-1;
                    Faux=np.flipud(F[0:imax2])
                    iumb=thresholdcross(Faux,umbral)
                    iumb=imax2-iumb+1
    
        if Pex==1:
            P1=iumb
    

    #% ---- P wave position ----
        if Pex==1:
            Faux=F[imax-1:QRS1]
            icero1=zerocross(Faux)
            icero1=imax+icero1-1
            Faux=np.flipud(F[0:imin])
            icero2=zerocross(Faux)
            icero2=imin-icero2+1
            Pp=round((icero1+icero2)/2) 
   
   #% ---- Check noise level ----
            inic=P1-40+1
            fin=P1-5
            if inic<=0:
                inic=1
            if fin<=0:
                fin=1
            Xaux=X[inic-1:fin]
            ruido=noiselevel(Xaux);
            if abs(Xpb[P1-1]-Xpb[Pp-1])<1.5*ruido and (Pp-P1)/Fs<40e-3:
                Pex=0
                type=2 
                P1=[]
                Pp=[]

    #% ---- P wave offset ----
        if Pex==1:
            umbral=(ymin)/Kpe
            Faux=F[imin-1:np.size(F)]
            iumb=thresholdcross(Faux,umbral)
            iumb=imin+iumb-1
            if iumb>=QRS1:
                ymin=min(F[imin-1:QRS1]) 
                iumb=min(F[imin-1:QRS1]) 
                iumb=imin+iumb-1
    
            P2=iumb
            if P2>=QRS1:
                P2=QRS1-1
    

#% ---- Check noise level ----
    if (Pex==1):
        Xaux=X[ibw-1:iew]
        ruido=noiselevel(Xaux)
        if abs(Xpb[Pp-1]-Xpb[P2-1])<=(1.5*ruido):
            Pex=0
            P1=[]
            Pp=[]
            P2=[]
            type=2
            

#% ---- Validation ----
    if (Pex==1):
        if P1>=P2 or Pp<=P1 or Pp>=P2 or P1<=prevt or (P2-P1)/Fs>180e-3 or (P2-P1)/Fs>150e-3:
            Pex=0
            P1=[]
            Pp=[]
            P2=[]
            type=2
        else:
            P1=[]
            Pp=[]
            P2=[]

    iew=iew-round(50e-3*Fs)
    ibw=ibw-round(50e-3*Fs)

    return P1,Pp,P2, type

def tbound (n,X,Xpa,F,PKni,Rp,QRS1,QRS2,PKnii,dermax,basel,RRm,Fs,ktb,kte,pco):
    
    T1=[]
    T2=[]
    Tp=[] 
    Tp2=[] 
    Ttype=6 
    peaks=0
    flag1=1
    kdis=round(50e-3*Fs)
    RRav=RRm/Fs
    M = []
    I = []
    M,I = peaksearch(F,0)
    if RRav>900e-3:
        itqlim=round(280e-3*Fs)
        ewindt=800e-3  
    elif RRav>800e-3:
        itqlim=round(250e-3*Fs)
        ewindt=600e-3   
    elif RRav>600e-3:
        itqlim=round(200e-3*Fs)
        ewindt=450e-3    
    else: 
        itqlim=round(170e-3*Fs)
        ewindt=450e-3   

    back=np.ones(1,6)
    bwindt=100e-3
    ibw=Rp+round(bwindt*Fs)
    iew=Rp+round(ewindt*Fs)

    if n>1:
        if RRav<750e-3:
            ibw=Rp+round(125e-3*Fs)
            iew=Rp+round(RRav*0.6*Fs)
        else:  
            ibw=Rp+round(bwindt*Fs)
            iew=Rp+round(ewindt*Fs)

    if ibw<=QRS2+kdis:
        iew=iew+QRS2-ibw+kdis
        ibw=QRS2+kdis

    if PKnii>0 and iew>PKnii-round(210e-3*Fs):
        iew=PKnii-round(210e-3*Fs)
    elif PKnii==0: 
        if Rp+round(400e-3*Fs)<=np.size(X):
            iew=Rp+round(400e-3*Fs)
        else:
            iew=np.size(X)
       

#% ---- Detecting upslope or downslope T ----
    ymax = max(F[ibw-1:iew])
    imax = max(F[ibw-1:iew])
    imax=imax+ibw-1
    ymin=min(F[ibw-1:iew])
    imin=min(F[ibw-1:iew])
    imin=imin+ibw-1
    if ymin>0 or ymax<0:
        if iew==PKnii-round(210e-3*Fs):
            if ymin>0: 
                ymin=0
            if ymax<0:
                ymax=0
        else:
            while (ymin>0 or ymax<0) and PKnii>0 and iew<PKnii-round(250e-3*Fs):
                iew=iew+round(25e-3*Fs)
                ymax=max(F[ibw-1:iew])
                imax=min(F[ibw-1:iew])
                imax=imax+ibw-1
                ymin=min(F[ibw-1:iew])
                imin=min(F[ibw-1:iew])
                imin=imin+ibw-1
        

    while (flag1==1)and(iew>ibw):
        peaks=0
        kint1=round(250e-3*Fs)
        kint2=round(300e-3*Fs)
        kend=50
        ampmi=0.075
        kk=3
        com=0.3
        if (-com*ymin<ymax and -ymin>com*ymax):
            if imin<imax:
                ymaxa=max(F[ibw-1:imin])
                imaxa=min(F[ibw-1:imin])
                imaxa=imaxa+ibw-1
                yminp=min(F[imax-1:iew])
                iminp=min(F[ibw-1:iew])
                iminp=iminp+imax-1
                if (ymaxa<ymax/kk and -yminp<-ymin/kk) or (ymaxa>=ymax/kk&-yminp>=-ymin/kk):
                    Ttype=1
                elif ymaxa>=ymax/kk or -yminp>=-ymin/kk:
                    peaks=1
            else:
                ymina=min(F[ibw-1:imax])
                imina=min(F[ibw-1:imax])
                imina=imina+ibw-1;
                ymaxp=max(F[imin-1:iew])
                imaxp=min(F[ibw-1:iew])
                imaxp=imaxp+imin-1
                if (ymaxp<ymax/kk&(-ymina)<(-ymin)/kk) or (ymaxp>=ymax/kk&(-ymina)>=(-ymin)/kk):
                    Ttype=0
                elif ymaxp>=ymax/kk or (-ymina)>=(-ymin)/kk:
                    peaks=1
        else: 
            peaks=1

#% ---- Different peaks ----
    if (peaks==1):
        if ymax>abs(ymin):
            ymina=min(F[ibw-1:imax])
            imina=min(F[ibw-1:imax])
            imina=imina+ibw-1
            yminp=min(F[imax-1:iew])
            iminp=min(F[ibw-1:iew])
            iminp=iminp+imax-1
            Faux=F[imina-1:np.size(F)]
            ncea=imina
            if ymina<0:
                ncea=zerocross(Faux)
                ncea=imina+ncea-1
       
            ampa=basel-X[ncea-1] 
            Faux=np.flipud(F[0:iminp])
            ncep=iminp
            if yminp<0:
                ncep=zerocross(Faux)
                ncep=iminp-ncep+1
      
            ampp=X[ncep-1]-basel
            if (ampa+ampp)>ampmi:
                if -ymina<ymax/pco and -yminp<ymax/pco:
                    Ttype=2
                elif -ymina>=ymax/pco and -yminp>=ymax/pco:
                    Ttype=4 
                elif -ymina>=ymax/pco and -yminp<ymax/pco:
                    Ttype=1 
                    ymin=ymina 
                    imin=imina
                elif -ymina<ymax/pco and -yminp>=ymax/pco:
                    Ttype=0
                    ymin=yminp
                    imin=iminp
         
            else: 
                Ttype=6
      
            
        elif ymax<abs(ymin):
            ymaxa=max(F[ibw-1:imin])
            imaxa = max(F[ibw-1:imin])
            imaxa=imaxa+ibw-1
            ymaxp=max(F[imin-1:iew])
            imaxp = max(F[imin-1:iew])
            imaxp=imaxp+imin-1
            Faux=F[imaxa-1:np.size(F)]
            ncea=imaxa;
            if ymaxa>0:
                ncea=zerocross(Faux)
                ncea=imaxa+ncea-1
      
            ampa=X[ncea-1]-basel
            Faux=np.flipud(F[0:imaxp])
            ncep=imaxp
            if ymaxp>0:
                ncep=zerocross(Faux)
                ncep=imaxp-ncep+1
      
            ampp=basel-X[ncep-1]
            if (ampa+ampp)>ampmi:
                if ymaxa<-ymin/pco and ymaxp<-ymin/pco:
                    Ttype=3 
                elif ymaxa>=-ymin/pco and ymaxp>=-ymin/pco:
                    Ttype=5
                elif ymaxa>=-ymin/pco & ymaxp<-ymin/pco:
                    Ttype=0 
                    ymax=ymaxa
                    imax=imaxa
                elif ymaxa<-ymin/pco and ymaxp>=-ymin/pco:
                    Ttype=1
                    ymax=ymaxp
                    imax=imaxp
         
            else:
                Ttype=6
       
      



#% ---- Normal or inverted T wave ----
    if (Ttype==0)or(Ttype==1):
        Tp2=[]
        if (Ttype==1): 
            yaux=ymax
            iaux=imax
            ymax=ymin
            imax=imin
            ymin=yaux
            imin=iaux
   
        umbral=(ymax)/ktb
        Faux=np.flipud(F[0:imax])
        iumba=thresholdcross(Faux,umbral)
        iumba=imax-iumba+1
    
        if (iumba<=QRS2):
            It=np.where(I<imax)[0]
            if not not It: 
                iumba=I[It[np.size(It)-1]]
                if iumba<=QRS2:
                    iumba=QRS2+2
      	
      
   
        T1=iumba
        if abs(ymin)>=0.41: 
            kte1=kte*2
        elif abs(ymin)>=0.35:
            kte1=kte*2-1
        elif abs(ymin)>=0.25: 
            kte1=kte*2-2
        elif abs(ymin)>=0.10: 
            kte1=kte*2-3
        elif abs(ymin)<0.10: 
            kte1=kte
   
        if kte/back[0]>=1.1:
            umbral=ymin*back[0]/kte1
        else:
            umbral=ymin/1.1
   
        Faux=F[imin-1:np.size(F)]
        iumbp=thresholdcross(Faux,umbral)
        iumbp=imin+iumbp-1
        T2=iumbp 
        Faux=np.flipud(F[0:imin])
        icero1=zerocross(Faux)
        icero1=imin-icero1+1
        Faux=F[imax-1:np.size(F)]
        icero2=zerocross(Faux)
        icero2=imax+icero2-1
        icero=round((icero1+icero2)/2)
        if (icero>=T2 or icero<=T1):
            icero=T1+round((T2-T1)/2)
   
        Tp=icero
        back[0]=back[0]*1.8
   
#% ---- Upslope or downslope T wave ----
    elif (Ttype==2)or(Ttype==3):
        T1=[]
        Tp2=[]
        if (Ttype==3):
            ymax=ymin
            imax=imin
        if abs(ymax)>=0.41:
            kte1=kte*2      
        elif abs(ymax)>=0.30:
            kte1=kte*2-1
        elif abs(ymax)>=0.20:
            kte1=kte*2-2
        elif abs(ymax)>0.10:
            kte1=kte*2-3
        elif abs(ymax)<=0.10:
            kte1=kte
   
        if kte/back(3)>=1.1:
            umbral=(ymax)*back[2]/kte1
        else:
            umbral=ymax/1.1
   
        Faux=F[imax-1:np.size(F)]
        iumbp=thresholdcross(Faux,umbral)
        iumbp=imax+iumbp-1
        T2=iumbp
        Faux=np.flipud(F[0:imax])
        icero=zerocross(Faux)
        icero=imax-icero+1
        It=np.where(I<(imax-kdis))[0]
        if not not It:
            ipic=I[It[np.size(It)-1]]
            Tp=max(ipic,icero)
            if Tp<=QRS2:
                Tp=QRS2+1
   
        back[2]=back[2]*1.8
 
    #% ---- Biphasic T wave ----
    elif (Ttype==4)or(Ttype==5):
        if (Ttype==5):
            ymina=ymaxa
            imina=imaxa
            ymax=ymin
            imax=imin
            yminp=ymaxp
            iminp=imaxp
   
        umbral=(ymina)/ktb
        Faux=np.flipud(F[0:imina])
        iumba=thresholdcross(Faux,umbral)
        iumba=imina-iumba+1
        if (iumba<=QRS2):
            It=np.where(I<imina)[0]
            if not not It:
                iumba=I[It[np.size(It)-1]]
                if iumba<=QRS2:
                    iumba=QRS2+2
        T1=iumba
        if abs(yminp)>=0.41: 
            kte1=kte*2    
        elif abs(yminp)>=0.30: 
            kte1=kte*2-1
        elif abs(yminp)>=0.20:
            kte1=kte*2-2
        elif abs(yminp)>=0.10:
            kte1=kte*2-3
        elif abs(yminp)<0.10:
            kte1=kte
   
        if kte/back(5)>=1.1:
            umbral=(yminp)*back[4]/kte1
        else:
            umbral=yminp/1.1
   
        Faux=F[iminp-1:np.size(F)]
        iumbp=thresholdcross(Faux,umbral)
        iumbp=iminp+iumbp-1
        T2=iumbp
        Faux=np.flipud(F[0:iminp])
        icero1=zerocross(Faux)
        icero1=iminp-icero1+1
        Tp=icero1
        Faux=F[imina-1:np.size(F)]
        icero2=zerocross(Faux)
        icero2=imina+icero2-1
        Tp2=icero2
        if Tp<Tp2: 
            Tp=Tp2
        back[4]=back[4]*1.8

#% ---- Validation ----
    if (T2-QRS1)<950e-3*Fs and((PKnii-T2>itqlim) or (PKnii==0) or (T2-QRS1<400e-3*Fs)):
        flag1=0
    else: 
        if iew>ibw+100e-3*Fs:
            iew=iew-round(50e-3*Fs)
        else:
            iew=iew-round(25e-3*Fs)
   
        ymin=min(F[ibw-1:iew]) 
        imin=min(F[ibw-1:iew]) 
        imin=ibw+imin-1
        ymax=max(F[ibw-1:iew])
        imax=max(F[ibw-1:iew])
        imax=ibw+imax-1


    if (PKnii-round(100e-3*Fs))<T2 and (PKnii!=0):
        Ttype=6

    if Ttype==6:
        T1=[]
        Tp2=[]
        Tp=[]
        T2=[]

    return T1,Tp2,Tp,T2,Ttype

def Qwave(n,X,D,Der,PKni,Rp,M,I,ymax,imax,ymin,imin,dermax,type,Sgran,Fs,Kq,Kr):
    
    Qp=[]
    QRS1=[]
    Qex=1
    crece=0
    Kq=1.8
    
    Daux=np.flipud(Der[0:imax+1])
    ncero=zerocross(Daux)
    ncero=imax+1-ncero+1
#   %if abs(Der(ncero+1))<abs(Der(ncero)) ncero=ncero+1; end
    Daux=np.flipud(D[0:imax])
    nceau=zerocross(Daux)
    nceau=imax-nceau+1
    #%if abs(D(nceau+1))<abs(D(nceau)) nceau=nceau+1; end 
    if not(nceau) or not(ncero):
        Qex=0
        
    if (Rp-ncero)/Fs>80e-3 and (Rp<=PKni):
        Qex=0

    if (Qex==1):
        Iq=np.where(I<nceau)[0]
        if not(Iq):
            Qex=0

    if Qex==1:
        mpic=I[Iq[np.size(Iq)-1]]

#% ---- Protection against cases in which the derivative almost crosses zero
#% ----
        Iq=np.where(I<imax)[0]
        icep=I[Iq[np.size(Iq)-1]]
        if  abs(D[mpic-1])>dermax/12 and not(icep>mpic and abs(D[icep-1])<dermax/50): 
       
 #    % ---- Detection of P wave joint to Q wave ----
            if (Rp-mpic)/Fs>90e-3 or ((nceau-mpic)/Fs>30e-3&Rp<=PKni):
                Qex=0        
       

     #% ---- Q wave onset ----
            umbral=D[mpic]/Kq;
            Daux=np.flipud(D[0:mpic])
            iumb=thresholdcross(Daux,umbral)
            iumb=mpic-iumb+1
            Iq=np.where(I<mpic)[0] 
            if not not(Iq):
                ipic=I[Iq[np.size(Iq)-1]]
                if ipic>iumb:
                    iumb=ipic
            
       
            if (Rp-iumb)/Fs>120e-3:
                Qex=0 
     

#% ---- Check that Q wave does not start with a small upslope peak ----
            ilimp=iumb-round(30e-3*Fs);
            if ilimp<=0:
                ilimp=1
            Daux=D[ilimp-1:iumb]
            ymin2=min(Daux)
            imin2=min(Daux)
            if not not(imin2):
                imin2=ilimp+imin2-1
            ymax2=max(Daux)
            imax2=max(Daux)
            if not not(imax2):
                imax2=ilimp+imax2-1
            if abs(ymin2)>=dermax/20:
                umbral=(ymin2)/Kq  
                Daux=np.flipud(D[0:imin2])
                iumb2=thresholdcross(Daux,umbral)
                iumb2=imin2-iumb2+1
                if (iumb2>iumb-round(40e-3*Fs)):
                    iumb=iumb2 
          
       

            if (imax2<imin2 or abs(ymin2)<dermax/20) and abs(ymax2)>dermax/20:
                umbral=(ymax2)/Kq 
                Daux=np.flipud(D[0:imax2])
                iumb2=thresholdcross(Daux,umbral)
                iumb2=imax2-iumb2+1
                if (iumb2>iumb-round(40e-3*Fs)):
                    iumb=iumb2
           
        
        else:
            crece=1     
   

#% ---- To detect high frequency Q waves, use the unfiltered derivative ----
        if (crece==1):
            Kq=Kq+1
            Md,Id=peaksearch(Der,0)
            Iq=np.where(Id<ncero)
            if not not(Iq):
                mpic=Id[Iq[np.size(Iq)-1]]
                mpic=testpeak(Der,mpic,Fs,0)
      
            if not(Iq) or abs(Der(mpic))<dermax/10 or (ncero-mpic)/Fs>30e-3:
                Qex=0   
      
            if (Qex==1):
                umbral=(Der(mpic))/2.8 
                Daux=np.flipud(Der[0:mpic])
                iumb=thresholdcross(Daux,umbral)
                iumb=mpic-iumb+1
                if (Rp-iumb)/Fs>80e-3:
                    Qex=0 
        

#% ---- If there is not Q wave, search for the onset of R wave ----
    if (Qex==0):
        if D[imax-1]>=4:
            Kr=Kr*6-2
        elif D(imax)>=3:
            Kr=Kr*4-2
        elif D(imax)>=1.5:
            Kr=Kr*2-2
    
        umbral=D(imax)/Kr
   
        Daux=np.flipud(D[0:imax])
        iumb=thresholdcross(Daux,umbral)
        if not(iumb):
            iumb=imax-round(30e-3*Fs)
        else: 
            iumb=imax-iumb+1
         

#% ---- Check that R wave does not start with a small peak ----
        ib=max(1,iumb-round(30e-3*Fs))
        Daux=D[ib-1:iumb]
        ymax2=max(Daux)
        imax2=max(Daux)
        if not not(imax2):
            imax2=ib+imax2-1
        if ymax2>=dermax/50:
            imax=imax2
            ymax=ymax2 
            umbral=ymax2/1.5
            Daux=np.flipud(D[0:imax])
            iumb2=thresholdcross(Daux,umbral)
            iumb2=imax-iumb2+1
            if iumb2>iumb-round(36e-3*Fs):
                iumb=iumb2 
       
    if (Qex==1):
        Qp=ncero 
    else:
        Qp=[]
    QRS1=iumb

    if not not(Qp):
        Qp=testpeak(X,Qp,Fs,0)
        
    return Qp, QRS1, type

def Rwave(n,X,Xpb,D,Der,PKni,M,I,Fs,Kr,Ks,Krr):

    QRS1=[]
    QRS2=[]
    Qp=[]
    Rp=[]
    Sp=[]
    R2p=[]
    Rex=0
    Qex=0
    Sex=0
    R2ex=0
    
    type=0
    noR=0
    Sgran=0 

#% ---- Previous and former peaks ----

    Ir=np.where(I>PKni)[0]
    mpicd=I[Ir[0]]
    Ir=np.where(I<PKni)[0]
    mpici=I[Ir[np.size(Ir)-1]]

#% ---- RSR' type? ----
    ydi=D[mpici]
    ydd=D[mpicd]
    ymaxaux=max(abs(ydd),abs(ydi))
    kpi=2

    if (Xpb[PKni-1]<0) or (ydi>0&ydd<0 and (kpi*ydi<(-1)*ydd or kpi*(-1)*ydd<ydi)):
        perc=0.25
        if (Xpb[PKni-1]>0 and ydi>0&ydd<0) or ((1+perc)*(-1)*ydi>ydd and(1-perc)*(-1)*ydi<ydd):
#      % ---- RSR'type ----
            type=2
            perc=0.35
            if (Xpb[PKni-1]<0):
#         % ---- PKni corresponds to S wave, R' will be on the right and R
#         % on the left ----
                Daux=D[mpicd-1:np.size(D)]
                ncero=zerocross(Daux)
                if not(ncero):
                    return
                ncero=mpicd+ncero-1
                Ir=np.where(I>ncero)[0]
                if not not(Ir):
                    mpda=I[Ir[0]]
                    if ((-1*D[mpda]<ydd/5) and (abs(Xpb[ncero])<abs(Xpb[PKni])/10)):  #% JGM
                         type=3 #% ---- Very large Q or S wave ----
                    else:
                        R2p=ncero
                        Sp=PKni
                        Daux=np.flipud(D[0:mpici])
                        ncero=zerocross(Daux)
                        if not(ncero): 
                            return
                        ncero=mpici-ncero+1
                        Ir=np.where(I<ncero)[0]
                        if not not(Ir):
                            mpda=I[Ir[np.size(Ir)-1]]
                            if ((D[mpda-1]<-ydi/5) and (abs(Xpb[ncero-1])<abs(Xpb[PKni-1])/10)): # % JGM
                                 type=3
                             # % ---- Very large Q or S wave ----
                            else:
                                 Rp=ncero
                             
            elif abs(ydi)<abs(ydd):
#         % ---- PKni corresponds to R wave ----
                 Daux=D[mpicd-1:np.size(D)]
                 ncero=zerocross(Daux)
                 if not(ncero):
                     return 
                 ncero=mpicd+ncero-1
                 Ir=np.where(I>ncero)[0] 
                 if not not(Ir): 
                     mpic=I[Ir[0]]
                     if not((1+perc)*abs(ydd)>abs(D(mpic)) and (1-perc)*abs(ydd)<abs(D[mpic-1])):
                         type=1# % ---- Normal QRS type ----
                     else: 
                         Sp=ncero
                         Daux=D[Sp-1:np.size(D)]
                         ncero=zerocross(Daux)
                         if not(ncero): 
                             return
                         ncero=Sp+ncero-1
                         R2p=ncero
                         Rp=PKni
            
            elif abs(ydi)>abs(ydd):
#            % ---- PKnii corresponds to R' wave ----
                Daux=np.flipud(D[0:mpici])
                ncero=zerocross(Daux)
                ncero=mpici-ncero+1
                if not not(ncero): 
                    return
                Ir=np.where(I<ncero)[0]
                if not not(Ir):
                    mpic=I[Ir[np.size(Ir)-1]]
                    if (not((1+perc)*abs(ydi)>abs(D[mpic-1])and(1-perc)*abs(ydi)<abs(D[mpic-1]))):
                        type=1# % ---- Normal QRS type ----
                    else: 
                        Sp=ncero
                        Daux=np.flipud(D[0:Sp])
                        ncero=zerocross(Daux)
                        if not(ncero):
                            return
                        ncero=Sp-ncero+1
                        Rp=ncero
                        R2p=PKni
                 
        if (type==2)and(R2p-Rp)/Fs>150e-3:
            if Xpb(PKni)>0:
                type=1#  % ---- Normal QRS type ----
            else:
                type=3# % ---- Very large Q or S wave ----
            
      
        else:
            type=3 #% ---- Very large Q or S wave ----  
   
    else:
        type=1#   % ---- Normal QRS type ----

   

#% ---- Onset and offset of RSR' ----
    if (type==2):
        Ir=np.where(I>R2p)[0]
        mpicd=I[Ir[0]]
        Ir=np.where(I<Rp)[0]
        mpici=I[Ir[np.size(Ir)-1]]
        R2p=testpeak(X,R2p,Fs,1)
        Sp=testpeak(X,Sp,Fs,0)
        Rp=testpeak(X,Rp,Fs,1)
        umbral=D[mpicd-1]/Krr
        Daux=D[mpicd-1:np.size(D-1)]
        QRS2=thresholdcross(Daux,umbral)
        QRS2=mpicd+QRS2-1
        if not(QRS2): 
            return
        umbral=X[QRS2-1]
        Xaux=np.flipud(X[0:R2p])
        S1=thresholdcross(Xaux,umbral)
        S1=R2p-S1+1
        umbral=D[mpici-1]/Kr
        Daux=np.flipud(D[0:mpici])
        QRS1=thresholdcross(Daux,umbral)
        QRS1=mpici-QRS1+1
        if not(QRS1):
            return
        umbral=X[QRS1-1]
        Xaux=X[Rp-1:np.size(X)]
        Q2=thresholdcross(Xaux,umbral)
        Q2=Rp+Q2-1

#% ---- R wave location in the very large Q or S wave case ----
    if (type==3):
        R2p=[]
        Rp=[]
        Sp=[]
        Daux=D[mpicd-1:np.size(D)]
        nrted=zerocross(Daux)
        nrted=mpicd+nrted-1
        Daux=np.flipud(D[0:mpici])
        nrtei=zerocross(Daux)
        nrtei=mpici-nrtei+1
        prr=1.4
        if (abs(D[mpicd-1])>prr*abs(D(mpici))or(PKni-nrtei)>(nrted-PKni)):
#         % ---- PKni corresponds to Q wave, R wave will be on the right
# % ----
             Daux=D[mpicd-1:np.size(D)]
             ncero=zerocross(Daux)
             if not(ncero):
                 return
             ncero=mpicd+ncero-1
             Ir=np.where(I>ncero)[0]
             mpda=I[Ir[0]]
             if (ncero-PKni)/Fs>150e-3 or (-1)*D[mpda-1]<ydd/10:
                 Sgran=1
             else:
                 Rp=ncero
         
        else:
            Sgran=1
      
        if Sgran==1:
#         % ---- PKni corresponds to S wave, R wave will be on the left ---- 
             Daux=np.flipud(D[0:mpici])
             ncero=zerocross(Daux)
             if not(ncero):
                 return
             ncero=mpici-ncero+1
             ilim=ncero-round(60e-3*Fs)
             if ilim<=0:
                 ilim=1
             Daux=D[ilim-1:ncero]
             if (not not(Daux)):
                 ymax2=max(Daux)
                 imax2=max(Daux)
                 imax2=ilim+imax2-1
       		  #%if (PKni-ncero)/Fs>140e-3|(ymax2<(-1)*ydi/10)|((Xpb(ncero))<-abs(Xpb(PKni))*1/6)
                 if (PKni-ncero)/Fs>140e-3 or (ymax2<(-1)*ydi/100) or (Xpb[ncero-1] < -abs(Xpb[PKni-1])*1/3):  	#%JGM
                     Rp=PKni  
                     type=4 # % ---- QS type ----
                 else:
                     Rp=ncero
#% ---- R wave location in the normal QRS type ----
    if (type==1):
        Rp=PKni
        R2p=[]

    return QRS1, Rp,Sp,R2p,QRS2,ymaxaux, type,Sgran

def Swave(n,X,D,Der,PKni,Rp,Sp,M,I,ymax,imax,ymin,imin,dermax,type,Sgran,Fs,Kr,Ks):
#% ---- S wave and QRS offset ----

    Sp=[]
    QRS2=[]
    Sex=1
    crece=0
    iumb=[]
    Daux=Der[imin-1:np.size(Der)]
    ncero=zerocross(Daux)
    ncero=imin+ncero-1
    #%if abs(Der(ncero-1))< abs(Der(ncero)) ncero=ncero-1; end
    Daux=D[imin-1:np.size(D)]
    nceau=zerocross(Daux)
    nceau=imin+nceau-1 
    #%if abs(D(nceau-1))<abs(D(nceau)) nceau=nceau-1; end

    if not(nceau) or not(ncero):
        Sex=0

    if (ncero-Rp)/Fs>130e-3&Rp>=PKni:
        Sex=0

    if nceau<PKni&(X[PKni-1]<0):
        ncero=PKni
        nceau=PKni

    if (Sex==1):
        if not not(Sp) and Sp==PKni:
            ilim=nceau+round(140e-3*Fs) 
        else:
            ilim=nceau+round(80e-3*Fs)

        if ilim>=np.size(D):
            ilim=np.size(D)
        Daux=D[nceau-1:ilim]
        ypic=max(Daux)
        mpic=max(Daux)
        mpic=nceau+mpic-1
        if ypic<dermax/10:
            Iq=np.where(I>=nceau)[0] 
            mpic=I[Iq[0]]
        Iq=np.where(I>imin)[0]
        icep=I[Iq[0]]

#% ---- Protection against cases in which the derivative almost exceeds zero
#% ----
        if abs(D[mpic-1])>dermax/30 and (not(icep<mpic and abs(D[icep-1])<dermax/50) or PKni==ncero): #%Antes /30.
            if (D[mpic-1])>=6.2:
                Ks=3*Ks+1
            elif (D[mpic-1])>=4.75:
                Ks=3*Ks
            elif (D[mpic-1])>=4:
                Ks=3*Ks-1 
      
            umbral=(D[mpic-1])/Ks
            inicio=mpic+round(10e-3*Fs)
            Daux=D[mpic-1:np.size(D)]
            iumb=thresholdcross(Daux,umbral)
            iumb=mpic+iumb-1
            Iq=np.where(I>inicio)[0]
            if not not(Iq):#                  %RBL
                ipic=I[Iq[0]]
                if (ipic<iumb) and D[ipic-1]<dermax/15:
                    iumb=ipic 
     
            if (iumb-Rp)/Fs>200e-3:
          
#        ---- There is not S wave ----
              umbral=(D[mpic-1])/Kr
              Daux=D[mpic-1:np.size(D)]
              iumb=thresholdcross(Daux,umbral)
              iumb=mpic+iumb-1 
              inicio=mpic+round(10e-3*Fs)
              Is=np.where(I>inicio)
              ipic=I[Is[0]]
              if ipic<iumb and D[ipic-1]<dermax/3:
                  iumb=ipic 
          
          
        else:        
#   % ---- Working with the unfiltered derivative ----
            Daux=Der[ncero-1:np.size(Der)]
            Md, Id=peaksearch(Der,0)
            Is=np.where(Id>ncero)[0]
            mpic=Id[Is[0]]
            mpic=testpeak(Der,mpic,Fs,1)
            if abs(Der[mpic-1])<dermax/10 and Rp>=PKni: 
                Sex=0 
    
            if (Sex==1):
                umbral=(Der[mpic-1])/Ks
                Daux=Der[mpic-1:np.size(Der)]
                iumb=thresholdcross(Daux,umbral)
                iumb=mpic+iumb-1 
                inicio=mpic+round(10e-3*Fs)
                Is=np.where(I>inicio)[0]
                ipic=I[Is[0]]
                if ipic<iumb:
                    iumb=ipic 
         
                if (iumb-Rp)/Fs>200e-3:
             
#       % ---- There is not S wave ----
                    umbral=(D[mpic-1])/Kr
                    Daux=D[mpic-1:np.size(D)]
                    iumb=thresholdcross(Daux,umbral)
                    iumb=mpic+iumb-1 
                    inicio=mpic+round(10e-3*Fs)
                    Is=np.where(I>inicio)[0]
                    ipic=I[Is[0]]
                    if ipic<iumb and D[ipic-1]<dermax/3:
                        iumb=ipic 
             
    if not not (iumb) and (iumb-Rp)/Fs>200e-3 and Rp>=PKni:
        Sex=0
        Sp=[]     
        #% ---- If there is not S wave, search for the onset of R wave ---- 
    if (Sex==0): 
        umbral=(D[imin-1])/Kr
        Daux=D[imin-1:np.size(D)]
        iumb=thresholdcross(Daux,umbral)
        iumb=imin+iumb-1
        inicio=imin+round(10e-3*Fs)
        Is=np.where(I>inicio)[0]
        ipic=I[Is[0]]
        if ipic<iumb:
            iumb=ipic 

    if (Sex==1):
        Sp=ncero
    else:
        Sp=[]
        
    QRS2=iumb
    if QRS2<PKni:
        QRS2=PKni+1
    if Sex==1:
        Sp=testpeak(X,Sp,Fs,0)   

    return Sp, QRS2, type, Sgran

#function [QRS1,Qp,Rp,Sp,R2p,QRS2,dermax,type,Sgran]=qrsbound(n,X,Xpb,D,Der,PKni,prevt,Fs,Kq,Kr,Ks,Krr)
def qrsbound(n,X,Xpb,D,Der,PKni,prevt,Fs,Kq,Kr,Ks,Krr):
    #% ---- QRS complex peak positions and limits depending on morphology ---- 
    #% ---- Initialization ----
    QRS1=[]
    QRS2=[]
    Qp=[]
    Rp=[]
    Sp=[]
    R2p=[]
    Rex=0
    Qex=0
    Sex=0
    R2ex=0
    type=0
    noR=0
    Sgran=0

    #% ---- R wave identification ----
    M,I=peaksearch(D,0)
    QRS1,Rp,Sp,R2p,QRS2,ymaxaux,type,Sgran=Rwave(n,X,Xpb,D,Der,PKni,M,I,Fs,Kr,Ks,Krr)

    Ir=np.where(I<Rp)[0]
    imax=I[Ir[np.size(Ir)-1]]
    ymax=M[Ir[np.size(Ir)-1]]
    Ir=np.where(I>Rp)[0] 
    imin=I[Ir[0]]
    ymin=M[Ir[0]]
    if type==2:
        dermax=max(abs(ymax),abs(ymin))
        #% ---- QRS type ----
        #% Protection against cases in which the derivative presents several
        #% peaks.
    if (type==1) or (type==3):
        if (ymax>ymaxaux/4):
            inicio=Rp-round(70e-3*Fs)
            Daux=D[inicio-1:Rp]
            ymaxa=max(Daux)
            imaxa=max(Daux)
            imaxa=inicio+imaxa-1
            ilim=Rp+round(70e-3*Fs)
            Daux=D[Rp-1:ilim]
            ymina=min(Daux)
            imina=min(Daux)
            imina=Rp+imina-1
            if ymaxa>ymax:
                ymax=ymaxa 
                imax=imaxa 
            if ymina<ymin:
                ymin=ymina
                imin=imina 

        dermax=max(abs(ymax),abs(ymin))
        ilim=imax-round(70e-3*Fs)
        ilim2=imax-round(30e-3*Fs)
        if ymax>ymaxaux/4:
            Daux=D[ilim-1:ilim2]
            ymaxa=max(Daux)
            imaxa=max(Daux)
            imaxa=ilim+imaxa-1
            if ymaxa>dermax/5:
                ymax=ymaxa
                imax=imaxa 
      
            ilim=imin+round(40e-3*Fs)
            ilim2=imin+round(100e-3*Fs)
            Daux=D[ilim-1:ilim2]
            ymina=min(Daux)
            ymina=min(Daux)
            imina=ilim+imina-1
            if abs(ymina)>dermax/5:
                ymin=ymina
                imin=imina 
   
    #% ---- QS type ----   
    if (type==4):  
        inicia=Rp-round(150e-3*Fs)
        Daux=D[inicia-1:Rp]
        ymin=min(Daux)
        imin=min(Daux)
        imin=inicia+imin-1
        ilim=Rp+round(180e-3*Fs) 
        Daux=D[Rp-1:ilim]
        ymax=max(Daux)
        imax=max(Daux)
        imax=Rp+imax-1
        dermax=max(abs(ymax),abs(ymin))

        umbral=ymin/Kr
        Daux=np.flipud(D[0:imin])
        QRS1=thresholdcross(Daux,umbral)
        QRS1=imin-QRS1+1

        ilim=QRS1-round(35e-3*Fs)
        Daux=D[ilim-1:QRS1]
        ymax2=max(Daux)
        imax2=max(Daux)
        imax2=ilim+imax2-1
        yaux=min(Daux)
        iaux=min(Daux)
        iaux=ilim+iaux-1

        if (ymax2)>=(dermax/30):
            Daux=np.flipud(D[0:imax2])
            umbral=ymax2/2
            iumb2=thresholdcross(Daux,umbral)
            iumb2=imax2-iumb2+1
            if iumb2>=QRS1-round(30e-3*Fs): 
                QRS1=iumb2
        
    
        if (abs(yaux)>=dermax/30) and (iaux<imax2):
            Daux=np.flipud(D[0:iaux])
            umbral=yaux/2
            iumb2=thresholdcross(Daux,umbral)
            iumb2=iaux-iumb2+1
            if iumb2>QRS1-round(50e-3*Fs): 
                QRS1=iumb2
        
        umbral=ymax/Kr
        Daux=D[imax-1:np.size(D)]
        QRS2=thresholdcross(Daux,umbral)
        QRS2=imax+QRS2-1
        ilim=Rp+round(180e-3*Fs)
        if (QRS2-QRS1)/Fs<80e-3:
            Daux=D[Rp-1:ilim]
            ymax2=max(Daux)
            imax2=max(Daux)
            imax2=Rp+imax2-1
            if ymax2>ymax:
                umbral=ymax2/Kr
                Daux=D[imax2-1:np.size(D)]
                QRS2=thresholdcross(Daux,umbral)
                QRS2=imax2+QRS2-1
        
        ilim=QRS2+round(20e-3*Fs)
        Daux=D[QRS2-1:ilim]
        ymin2=min(Daux) 
        imin2=min(Daux) 
        if not not(imin2):
            imin2=QRS2+imin2-1
            yaux=max(Daux)
            iaux=max(Daux)
            iaux=QRS2+iaux-1
            if abs(ymin2)>dermax/20:
                umbral=ymin2/2
                Daux=D[imin2-1:np.size(D)]
                iumb2=thresholdcross(Daux,umbral)
                iumb2=imin2+iumb2-1
                if iumb2<QRS2+round(30e-3*Fs):
                    QRS2=iumb2
          
    #% ---- QRS type ----
    if (type==1)or(type==3):
    #% ---- Q wave and QRS onset ----
        Qp, QRS1, type=Qwave(n,X,D,Der,PKni,Rp,M,I,ymax,imax,ymin,imin,dermax,type,Sgran,Fs,Kq,Kr)

    #% ---- S wave and QRS offset ----
        Sp,QRS2,type,Sgran=Swave(n,X,D,Der,PKni,Rp,Sp,M,I,ymax,imax,ymin,imin,dermax,type,Sgran,Fs,Kr,Ks)

    if QRS1<prevt:
        QRS1=prevt+2
        
    return QRS1, Qp, Rp, Sp, R2p, QRS2, dermax, type,Sgran

def proces(fid,X,Xpa,Xpb,D,F,Der,ti,tf,iprimerqrs,nqrs,iqrs,atyp,ns,Fs,nl,res,prewindt,Kq,Kr,Ks,Krr,Kpb,Kpe,Ktb,Kte,pco):
    
    Kr=5
    ns=np.size(X)
#    T=linspace(ti,tf,ns)

    #% ---- Initialization of RR interval ----
    if nqrs>=1:
        RRm=(iqrs[nqrs-1]-iqrs[0])/(nqrs-1)

    #% ---- Initializations ----
    ipbeg=np.zeros(1,nqrs)
    ippos=np.zeros(1,nqrs)
    ipend=np.zeros(1,nqrs)
    iqbeg=np.zeros(1,nqrs)
    iqpos=np.zeros(1,nqrs)
    iqend=np.zeros(1,nqrs)
    irpos=np.zeros(1,nqrs)
    ir2pos=np.zeros(1,nqrs)
    isbeg=np.zeros(1,nqrs)
    ispos=np.zeros(1,nqrs)
    isend=np.zeros(1,nqrs)
    itbeg=np.zeros(1,nqrs)
    itpos=np.zeros(1,nqrs)
    it2pos=np.zeros(1,nqrs)
    itend=np.zeros(1,nqrs)
    Qamp=np.zeros(1,nqrs)
    Ramp=np.zeros(1,nqrs)
    Samp=np.zeros(1,nqrs)
    iRR=np.diff(iqrs)
    
    POSPonset =[]
    POSP =[]
    POSPoffset = []
    POSQRSonset = []
    POSQ = []
    POSR = []
    POSfiducial = []
    POSS = []
    POSR2 = []
    POSQRSoffset = []
    POSTonset = []
    POST = []
    POST2 = []
    POSToffset = []

    AMPP = []
    AMPQ = []
    AMPR = []
    AMPS = []
    AMPR2 = []
    AMPT = []
    AMPT2 = []

    POSPonset[0:nqrs]=math.nan
    POSP[0:nqrs]=math.nan
    POSPoffset[0:nqrs]=math.nan
    POSQRSonset[0:nqrs]=math.nan
    POSQ[0:nqrs]=math.nan
    POSR[0:nqrs]=math.nan
    POSfiducial[0:nqrs]=math.nan
    POSS[0:nqrs]=math.nan
    POSR2[0:nqrs]=math.nan
    POSQRSoffset[0:nqrs]=math.nan
    POSTonset[0:nqrs]=math.nan
    POST[0:nqrs]=math.nan
    POST2[0:nqrs]=math.nan
    POSToffset[0:nqrs]=math.nan

    AMPP[0:nqrs]=math.nan
    AMPQ[0:nqrs]=math.nan
    AMPR[0:nqrs]=math.nan
    AMPS[0:nqrs]=math.nan
    AMPR2[0:nqrs]=math.nan
    AMPT[0:nqrs]=math.nan
    AMPT2[0:nqrs]=math.nan

    POS_QT=np.zeros(nqrs,1)
    VAL_QT=np.zeros(nqrs,1)
    VAL_QTC=np.zeros(nqrs,1)

    POS_QRS=np.zeros(nqrs,1)
    VAL_QRS=np.zeros(nqrs,1)

    POS_ANNOT=np.zeros(nqrs*10,1)
    ANNOT=np.zeros(nqrs*10,1)
    NUMFIELD=np.zeros(nqrs*10,1)
    SUBTYPEFIELD=np.zeros(nqrs*10,1)
    CHANFIELD=np.zeros(nqrs*10,1)

    n=1
    a=1
    q=1

    #% ---- Beat processing loop ----
    while n<nqrs-1 and (ns-iqrs[n-1])/Fs>500e-3:
        if n>=iprimerqrs and n<nqrs-1:
            basel=0
         
            #% ---- Analysis window ----
            bwind=max(iqrs[n-1]-round(RRm),iqrs[n-1]-round(Fs))
            if bwind<0:
                bwind=0
            if n<nqrs:
                ewind=max(iqrs[n-1]+round(RRm),iqrs[n])
            else:
                ewind=min(iqrs[n-1]+round(RRm),ns)
         
            if ewind>np.size(X):
                ewind=np.size(X)
         
            PKni=iqrs[n-1]-bwind
            if n<nqrs:
                PKnii=iqrs[n]-bwind 
            else:
                PKnii=0
         
            if n>=iprimerqrs+1:
                prevt=itend[n-2]-bwind
            else:
                prevt=prewindt-ti-bwind
         
            if prevt<=0:
                prevt=0
             
        #% ---- Detection of QRS complex position and limits ----
            QRS1,Qp,Rp,Sp,R2p,QRS2,dermax,Rtype,Sgran=qrsbound(n,X[bwind:ewind],Xpb[bwind:ewind],D[bwind:ewind],Der[bwind:ewind],PKni,prevt,Fs,Kq,Kr,Ks,Krr)
            if not not(QRS1):
                iqbeg[n-1]=bwind+QRS1
            if not not(Qp):
                iqpos[n-1]=bwind+Qp
            if not not(Rp):
                irpos[n-1]=bwind+Rp
            if not not(Sp):
                ispos[n-1]=bwind+Sp
            if not not(R2p):
                ir2pos[n-1]=bwind+R2p
            if not not(QRS2):
                isend[n-1]=bwind+QRS2

        #% ---- P wave detection ----
            P1,Pp,P2,Ptype=pbound(n,X[bwind:ewind],Xpb[bwind:ewind],F[bwind:ewind],PKni,Rp,QRS1,prevt,dermax,Fs,Kpb,Kpe)
            if not not(P1):
                ipbeg[n-1]=bwind+P1
            if not not(Pp):
                ippos[n-1]=bwind+Pp
            if not not(P2):
                ipend[n-1]=bwind+P2
                
        #% ---- Baseline estimation ----
            nqui=round(15e-3*Fs)
            ntre=round(30e-3*Fs)
            ntre_q=round(10e-3*Fs)
            nqui_q=round(5e-3*Fs)
            if not not(P2):
                if (QRS1-P2)/Fs>33e-3:
                    Xaux=X[ipend[n]+nqui-1:iqbeg[n]-nqui]
                    basel=sum(Xaux)/np.size(Xaux)
                elif (QRS1==P2):
                    basel=X[iqbeg[n-1]]
                else:
                    Xaux=X[ipend[n]-1:iqbeg[n]]
                    basel=sum(Xaux)/np.size(Xaux)
                  
             
            else:
                Xaux=X[iqbeg[n]-ntre_q-nqui_q-1:iqbeg[n]-nqui_q]
                basel=sum(Xaux)/np.size(Xaux)
         

    #% ---- Q and S wave offset ----
            if not not(Qp) and (Rtype==1 or Rtype==3):
                if X[irpos[n-1]]>0:
                    Xaux=X[iqpos[n]-1:isend[n]]
                    Q2=thresholdcross(Xaux,basel)
                    if not not(Q2):
                        iqend[n-1]=irpos[n-1] 
                    else:
                        iqend[n-1]=iqpos[n-1]+Q2-1 
               
                else:
                    iqend[n-1]=irpos[n-1]
            
        
            if not not(Sp) and (Rtype==1 or Rtype==3):
                if X[irpos[n-1]]>0:
                    Xaux=np.flipud(X[iqbeg[n]-1:ispos[n]])
                    S1=thresholdcross(Xaux,basel) 
                    if not(S1):
                        isbeg[n-1]=irpos[n-1]
                    else:
                        isbeg[n-1]=Sp-S1+1 
                
                else:
                    isbeg[n-1]=irpos[n]-1
             
        
        #% ---- Mean RR interval ----
            RRm=calcrr(RRm,n,Fs,iqrs)

#% ---- T wave location and limits ----
            T1,Tp2,Tp,T2,Ttype=tbound(n,X[bwind:ewind],Xpa[bwind:ewind],F[bwind:ewind],PKni,Rp,QRS1,QRS2,PKnii,dermax,basel,RRm,Fs,Ktb,Kte,pco) 

            if not not(T1):
                itbeg[n-1]=bwind+T1
            if not not(Tp2):
                it2pos[n-1]=bwind+Tp2
            if not not(Tp):
                itpos[n-1]=bwind+Tp
            if not not(T2):
                itend[n-1]=bwind+T2

    #% ---- Q, R, S and R' wave amplitudes ----
            if (irpos[n-1]!=0 and Rtype!=4):
                irpos[n-1]=testpeak(X,irpos[n-1],Fs,1)# %JGM
            elif (irpos(n)!=0 and Rtype==4):
                irpos[n-1]=testpeak(X,irpos[n-1],Fs,0)# %JGM
         
            if ispos[n-1]!=0: 
                ispos[n-1]=testpeak(X,ispos[n-1],Fs,0)# %JGM 
         
            if ir2pos[n-1]!=0: 
                ir2pos[n-1]=testpeak(X,ir2pos[n-1],Fs,1)# %JGM 
         
            if not not(Qp):
                Qamp[n-1]=X[iqpos[n-1]]-basel
                AMPQ[n-1]=Qamp[n-1]
         
            if not not(Rp):
                Ramp[n-1]=X[irpos[n-1]]-basel
                AMPR[n-1]=Ramp[n-1]
         
            if not not(Sp) and Rtype!=4:
                Samp[n-1]=X[ispos[n-1]]-basel
                AMPS[n-1]=Samp[n-1]
         
            if not not(R2p):
                R2amp = []
                R2amp[n-1]=X[ir2pos[n-1]]-basel
                AMPR2[n-1]=R2amp[n-1]

    #% ---- P and T wave amplitudes ----
            if not not(Pp):
                AMPP[n-1]=X[ippos[n-1]]-basel
            if not not(Tp2):
                AMPT2[n-1]=X[it2pos[n-1]]-basel
            if not not(Tp):
                AMPT[n-1]=X[itpos[n-1]]-basel
            
        n += 1
    
    #     I = []
    #     I=np.where(POS_ANNOT!=0)
    #     if not not(I):
    #         POSANNOT=POSANNOT[I-1]
    #         ANNOT=ANNOT[I-1]
    #         NUMFIELD=NUMFIELD(I);
    #     SUBTYPEFIELD=SUBTYPEFIELD(I);
    #     CHANFIELD=CHANFIELD(I);
    # else POS_ANNOT=[];
    #      ANNOT=[];
    #      NUMFIELD=[];
    #      SUBTYPEFIELD=[];
    #      CHANFIELD=[];
    # end
    
    
    # I=find(VAL_QT~=0);
    # if ~isempty(I)
    #     VAL_QT=VAL_QT(I);
    #     POS_QT=POS_QT(I);
    #     VAL_QTC=VAL_QTC(I);
    #     iRRc=iRRc(I);
    #     VAL_QT=(VAL_QT./Fs)*1000;  % Value in ms.
    #     POS_QT=POS_QT./Fs;
    #     VAL_QTC=VAL_QT./(sqrt(iRRc))';
    # else VAL_QT=[];
    #      POS_QT=[];
    #      VAL_QTC=[];
    # end
    
    # I=find(VAL_QRS~=0);
    # if ~isempty(I)
    #    VAL_QRS=VAL_QRS(I);
    #    POS_QRS=POS_QRS(I);
    #    POS_QRS=POS_QRS./Fs;
    #    VAL_QRS=(VAL_QRS./Fs)*1000; % Value in ms.
    # else VAL_QT=[];
    #      POS_QRS=[];
    # end
    
    # I=find(Qamp~=0);
    # if ~isempty(I)
    #     AMP_Q=Qamp(I)'; 
    #     POS_Q=((iqpos(I)+ti)./Fs)';
    # else AMP_Q=[]; POS_Q=[];
    # end
    
    # I=find(Ramp~=0);
    # if ~isempty(I)
    #    if Rtype~=4
    #     AMP_R=Ramp(I)';
    #     POS_R=((irpos(I)+ti)./Fs)';
    #     Is=find(Samp~=0);
    #     if ~isempty(Is)
    #         AMP_S=Samp(Is)';
    #         POS_S=((ispos(Is)+ti)./Fs)';
    #      else AMP_S=[]; POS_S=[];
    #      end
    #    else AMP_S=Ramp(I)';
    #     POS_S=((irpos(I)+ti)./Fs)';
    #     AMP_Q=[]; POS_Q=[];
    #     AMP_R=[]; POS_R=[]; 
    # end
    # else
    #     AMP_S=[]; POS_S=[];
    #     AMP_R=[]; POS_R=[];
    # end
    
    # I=[iprimerqrs:nqrs-2];
    # POS.Ponset=POS.Ponset(I);
    # POS.P=POS.P(I);
    # POS.Poffset=POS.Poffset(I);
    # POS.QRSonset=POS.QRSonset(I);
    # POS.Q=POS.Q(I);
    # POS.R=POS.R(I);
    # POS.fiducial=POS.fiducial(I);
    # POS.S=POS.S(I);
    # POS.R2=POS.R2(I);
    # POS.QRSoffset=POS.QRSoffset(I);
    # POS.Tonset=POS.Tonset(I);
    # POS.T=POS.T(I);
    # POS.T2=POS.T2(I);
    # POS.Toffset=POS.Toffset(I);
    
    # AMP.P=AMP.P(I);
    # AMP.Q=AMP.Q(I);
    # AMP.R=AMP.R(I);
    # AMP.S=AMP.S(I);
    # AMP.R2=AMP.R2(I);
    # AMP.T=AMP.T(I);
    # AMP.T2=AMP.T2(I);
