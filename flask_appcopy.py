from gammalkod.petprocessing import petcalc
import numpy as np
import gammalkod.Solweig_v2015_metdata_noload as metload
import gammalkod.clearnessindex_2013b as ci
from noahtingb7.petprocessing import petcalc as petcalc1
import noahtingb7.clearnessindex_2013b as ci1
import requests
import json
import base64
import pandas as pd
import noahtingb8.Processing as pp

class AnalyticBase1:
    def getdata(self,cor=[57,12],types=["temperature_2m"],forcastdays=14):
            return requests.get('https://api.open-meteo.com/v1/forecast?latitude='+str(cor[0])+'&longitude='+str(cor[1])+'&hourly'+"=,".join(types)+"&forecast_days="+str(forcastdays))
    def gethourlyparams1(self):
        return ["temperature_2m", "relative_humidity_2m", "cloud_cover", "cloud_cover_low", "cloud_cover_mid", "cloud_cover_high", "wind_speed_10m", "wind_gusts_10m", 
                "shortwave_radiation", "direct_radiation", "diffuse_radiation"]
    def gethourlyparams2(self):
        return ["apparent_temperature", "surface_pressure", "snowfall", "rain", "showers", "weather_code", "snow_depth", "visibility","global_tilted_irradiance","direct_normal_irradiance","wind_direction_10m","wind_direction_10m"]

    def dumpadatan(self,dat,filename="data//data.json"):
        with open(filename,"w") as file:
            json.dump(dat,file)
    def readdatan(self,filename): 
        with open(filename,"r") as file:
            filen=file.read()
        return json.loads(filen)

#Skapar en ny fil med väderdata baserat på koordinater och plastsnamn
    def createnewfilewithparams(self,cordinater=[12,58],place="Nowere",timea=[2025]):
        print(place,cordinater)
        datan=self.getdata(cordinater[0],self.AB.gethourlyparams1(),14)
        filenameis="data/"+place+"_"+"_".join(["_".join(i.split(":")) for i in timea])
        with open("tempimpfile/latestfile.json","r") as file:
            listan=json.loads(file.read())     
        with open("./tempimpfile/latestfile.json","w") as file:
            if listan.get(place)==None:
                listan[place]=[filenameis+".json"]
            else:
                listan[place].append(filenameis+".json")
            json.dump(listan,file)  
        self.dumpadatan(datan,filename=filenameis+".json")

class Analyticuppermid1:
#Denna klass använder Analyticclowermid för att hantera och visulaisera väderdata


#Initialiserar klassen och skapar en instans av Analyticclowermid
    def __init__(self,filen="./data/Bergstena_Wed_Jan_29_23_50_59_2025.json"):
        self.AB=AnalyticBase1()
        self.filen=filen
        self.data=self.AB.readdatan(filen)
    
#Byter till den senaste filen
    def changetolatestfile(self):#inactive
        with open("./tempimpfile/latestfile.json","r") as file:
            lista=json.loads(file.read())
        if type(lista)==type([]) and len(lista)>1:
            filenn=lista[len(lista)-2]
            if filenn[len(filenn)-4:]==".json":
                self.filen=filenn
                self.data=self.AB.readdatan(self.filen)
            else:
                pass

#Byter till en specifik fil
    def changefile(self,filen):
        self.filen=filen
        self.data=self.AB.readdatan(self.filen)



def prognose(method,form):
    if method == "GET":
        return 
    if method == "POST":
        city = form["city"]
        print("running")
        obj=AnalyticBase1()
        data=obj.getdata([57,12],obj.gethourlyparams1(),1)
        print(data)
        print(json.loads(data))
        return 
        timestamp=data[""]
        # putting data in separate vectors
        veclen = timestamp.__len__()
        year = np.empty(veclen, int)
        month = np.empty(veclen, int)
        day = np.empty(veclen, int)
        hour = np.empty(veclen, int)
        minu = np.empty(veclen, int)
        year = np.empty(veclen, int)
        Ta = np.empty(veclen, float)
        RH = np.empty(veclen, float)
        radD = np.empty(veclen, float)
        radI = np.empty(veclen, float)
        radG = np.empty(veclen, float)
        Ws = np.empty(veclen, float)

        for i in range(0, veclen):
            year[i] = int(timestamp[i][0:4])
            month[i] = int(timestamp[i][5:7])
            day[i] = int(timestamp[i][8:10])
            hour[i] = int(timestamp[i][11:13])
            minu[i] = int(timestamp[i][14:16])
            Ta[i] = float(dict_loaded['data_vars']['air_temperature']['data'][i])
            RH[i] = float(dict_loaded['data_vars']['relative_humidity']['data'][i])
            radD[i] = float(dict_loaded['data_vars']['downward_diffuse']['data'][i])
            radI[i] = float(dict_loaded['data_vars']['downward_direct']['data'][i])
            Ws[i] = np.sqrt(float(dict_loaded['data_vars']['eastward_wind']['data'][i])**2 + float(dict_loaded['data_vars']['northward_wind']['data'][i])**2)
            
        with np.errstate(invalid='ignore'):
          radI[radI < 0.] = 0.
          radD[radD < 0.] = 0.
        radG = radD + radI

        poi_save = petcalcprognose(Ta, RH, Ws, radG, radD, radI, year, month, day, hour, minu, lat, lon, UTC)

        tab = pd.DataFrame(poi_save[1:,[1,2,22,24,26,33]])
        tab.columns = ['Day of Year', 'Hour','T_air','RH','Tmrt', 'PET']

        tabhtml = tab.to_html(classes='data', header="true")

        doy = poi_save[1:, 1]
        hour = poi_save[1:, 2]
        #petres = str(round(poi_save[:,26], 1))

        return print("petprognoseresult.html", result1=doy, result2=hour, result3=tabhtml)
#prognose("POST",{"city":"Goteborg"})
def index(form):
        month = int(form["month"])
        day = int(form["day"])
        hour = int(form["hour"])
        year = int(form["year"])
        minu = 30
        Ta = float(form["Ta"])        
        RH = float(form["RH"])
        Ws = float(form["Ws"])
        location = form["loc"]

        if month > 12 or month < 0:
            print("petresult.html","Incorrect month filled in")
        if day > 31 or day < 0:
            print("petresult.html","Incorrect day filled in")
        if hour > 23 or hour < 0:
            print("petresult.html","Incorrect hour filled in")
        if Ta > 60 or Ta < -75:
            print("petresult.html", "Unreasonable air temperature filled in",Ta)
        if RH > 100 or RH < 0:
            print("petresult.html", "Unreasonable relative humidity filled in")
        if Ws > 100 or Ws < 0:
            print("petresult.html", "Unreasonable Wind speed filled in")

        #day of year
        if (year % 4) == 0:
            if (year % 100) == 0:
                if (year % 400) == 0:
                    leapyear = 1
                else:
                    leapyear = 0
            else:
                leapyear = 1
        else:
            leapyear = 0

        if leapyear == 1:
            dayspermonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        else:
            dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        doy = np.sum(dayspermonth[0:month - 1]) + day

        # Currently looked to Gothenburg
        UTC = 0

        # Radiation
        P = -999.
        radG = 40.
        #print("loc",location)
        metdata = np.zeros((1, 24)) - 999.
        metdata[0, 0] = year
        metdata[0, 1] = doy
        metdata[0, 2] = hour
        metdata[0, 3] = minu
        metdata[0, 11] = Ta
        metdata[0, 10] = RH
        #print(metdata)

        _, altitude, _, zen, jday, _, _, _ = metload.Solweig_2015a_metdata_noload(metdata, location, UTC)
        #print(altitude,zen,jday)
        if altitude > 0.:
            I0, _, _, _, _ = ci.clearnessindex_2013b(zen, jday, Ta, RH / 100., 40, location, P)
            radG = I0*1#(1-form["c1"]*1e-2)*(1-form["c2"]*1e-2)*(1-form["c3"]*1e-2)
        else:
            radG = 0.

        # Main calculation
        if Ta is not None and RH is not None and Ws is not None and radG is not None:
            #print(Ta, RH, Ws, radG, year, month, day, hour, minu)           

            Tmrt, resultPET, _ = petcalc(Ta, RH, Ws, radG, int(year), month, day, hour, minu,location)
            result = str(round(resultPET, 1))
        radG1=radG


        # Main calculation
        if Ta is not None and RH is not None and Ws is not None and radG is not None:
            #print(Ta, RH, Ws, radG, year, month, day, hour, minu)           
            
            Tmrt1, resultPET1, _ = petcalc1(Ta, RH, Ws, 0, year, month, day, hour, minu,location)
            result1 = str(round(resultPET1, 1))
        params={"year":year,"month":month,"day":day,"hour":hour,"Ta":Ta,"SwR":0.,"RH":RH,"Ws":Ws,"sky":"Clear (100%)","radI":0., "radD":0., "loc":location,"UTC":0}
        resultPET2,Tmrt2=pp.indexflask(params)
        #print("petresult", result,result1)
        return resultPET,resultPET1,resultPET2,Tmrt,Tmrt1,Tmrt2

"""
import sun_position as sp
import matplotlib.pyplot as plt

timmm=np.linspace(0,28*24*600,28*24*60)
results1=np.array([.0]*28*24*60)
results2=np.array([.0]*28*24*60)

for i in range(28*24*60):
    time={"year":2025,"month":1,"day":(i//(24*60))%31,"hour":(i//60)%24,"min":i%60,'sec':0,"UTC":1}
    location = {'longitude': 12.0, 'latitude': 57.7, 'altitude': 3.}
    res=sp.sun_position(time,location)
    results1[i],results2[i]=float(res['zenith'][0]),float(res['azimuth'][0])
    print(results2[i])
print(results1)
print(results2)
plt.plot(timmm,results1)
plt.plot(timmm,np.abs(180-results2))
plt.show()"
"""
#index("POST",{"year":2025,"month":1,"day":12,"hour":12,"Ta":20,"RH":50,"Ws":20,"sky":"Clear (100%)","c1":2,"c2":95,"c3":2})
#print("[[10.79787497]] [[1.38233786]] [[12.]]\n 20.0 50.0 20.0 [[9.92866941]] 2025 1 12 12 30\n petresult.html 13.1")
