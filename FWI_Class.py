import math

class FWI:
    def __init__(self,temp,rhum,wind,precp):
        self.t=temp
        self.h=rhum
        self.w=wind
        self.p=precp

    def FFMC(self,fo):
        mo=147.2 * (101.0-fo)/(59.5+fo)
        mint=42.5*rf*math.exp(-100.0/(251.0-mo))*(1.0-math.exp(-6.93/rf))

        if (self.p > 0.5):
            rf=self.p-0.5
            if(mo<=150.0):
                mo=mo+mint
            elif (mo>150):
                mo=mo+mint+0.0015*((mo -150.0)**2)*(rf**0.5)

        if mo>250.0:
            mo=250.0
        ed=0.942*(self.h**0.679)+11.0*math.exp((self.h-100.0)/10.0)+0.18*(21.1-self.t)*(1-math.exp(-0.115*self.h))

        if(mo>ed):
            ko=0.424*(1.0-(self.h/100.0)**1.7)+0.0694*(self.w**0.5)*(1.0-(self.h/100.0)**8)
            kd=ko*(0.581*math.exp(0.0365*self.t))
            m=ed+(mo-ed)*(10**(-kd))

        elif(mo<ed):
            ew=0.618*(self.h**0.753)+10.0*math.exp((self.h-100.0)/10.0)+0.18*(21.1-self.t)*(1-math.exp(-0.115*self.h))
            if(mo<ew):
                ki=0.424*(1.0-((100-self.h)/100)**1.7)+0.0694*(self.w**0.5)*(1.0-((100-self.h)/100)**8.0)
                kw=ki*0.581*math.exp(0.0365*self.t)
            elif(mo>ew):
                m=mo


        elif(mo==ed):
            m=mo

        ffmc_val=59.5*(250.0-m)/(147.2+m)

        if(ffmc_val<=0.0):
            ffmc_val=0.0

        if(ffmc_val>101.0):
            ffmc_val=101.0

        return ffmc_val
    def DMC(self,po,mth):
        le=[6.5,7.5,9.0,12.8,13.9,13.9,12.4,10.9,9.4,8.0,7.0,6.0]
        t=self.t
        if(t<-1.1):
            t=-1.1

        if(self.p>1.5):
            re=0.92*self.p-1.27
            mo=20.0+math.exp((5.6348-po)/43.43)
            if(po<=33.0):
                b=100.0/(0.5+0.3*po)
            elif(po>33.0):
                if(po<=65):
                    b=14.0-(1.3*math.log(po))
                elif(po>65):
                    b=6.2*math.log(po)-17.2
            mr=mo+1000.0*re/(48.77*b*re)
            po=244.72-43.43*math.log(mr-20.0)

        elif(self.p<=1.5):
            pr=po

        k=1.894*(t+1.1)*(100-self.h)*le[mth-1]*(10.0**-6)
        dmc_val=po+100*k
        if(dmc_val<=1.0):
            dmc_val=1.0
        return dmc_val

    def DC(self,do,mth):
        lf=[-1.6,-1.6,-1.6,0.9,3.8,5.8,6.4,5.0,2.4,0.4,-1.6,-1.6]
        t=self.t
        if(t<-2.8):
            t=-2.8

        if(self.p>2.8):
            rd=0.83*self.p-1.27
            qo=800.0*math.exp(do/400.0)
            qr=qo+3.937*rd
            dr=400.0*math.log(800.0/qr)
            v=0.36*(t+2.8)+lf[mth-1]
            if(v<0.0):
                v=0.0
            dc_val=dr+0.5*v
        elif(self.p<=2.8):
            dc_val=(0.36*(t+2.8)+lf[mth-1])+do
        return dc_val

    def ISI(self,ffmc):
        mo=147.2 * (101.0-ffmc)/(59.5+ffmc)
        fw=math.exp(0.005039*self.w)
        ff=91.9*math.exp(-0.1386)*((1+m**5.31)/(4.93*(10**7)))
        isi_val=0.208*ff*fw
        return isi_val

    def BUI(self,dmc,dc):
        if (dmc<=0.4*dc):
            bui_val=(0.8*dmc*dc)/(dmc+0.4*dc)
        elif(dmc>0.4):
            bui_val=dmc-((1.0-0.8*dc)/(dmc+0.4*dc))*(0.92+(0.0114*dmc)**1.7)
        else:
            bui_val=0.0
        return bui_val

    def fwi(self,isi,bui):
        if(bui<=80):
            fd=(0.626*bui)**0.809 + 2.0
        elif(bui>80):
            fd=1000.0/(25.0+108.64*math.exp(-0.023*bui))

        b=0.1*isi*fd
        if(b>1):
            fwi=2.72*((0.434*math.log(b))**0.647)
        elif(b<=1):
            fwi=b

        return fwi

    def dsr(self,fwi):
        dsr=0.0272*(fwi**1.77)
        return dsr
