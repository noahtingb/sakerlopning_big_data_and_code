import flask_appcopy as flaskapp
import requests
import json
import numpy as np
import time
import matplotlib.pyplot as plt

def getdata(cor=[57,12],types=["temperature_2m"],forcastdays=14):
        return requests.get('https://api.open-meteo.com/v1/forecast?latitude='+str(cor[0])+'&longitude='+str(cor[1])+'&hourly='+",".join(types)+"&forecast_days="+str(forcastdays))
def gethourlyparams1():        
    return ["temperature_2m", "relative_humidity_2m", "cloud_cover", "cloud_cover_low", "cloud_cover_mid", "cloud_cover_high", "wind_speed_10m", "wind_gusts_10m", "shortwave_radiation", "direct_radiation", "diffuse_radiation"]

def filehandeling():
    try:
        f = json.loads("filename.json")
    except FileNotFoundError:
        result=getdata(cor=[loca["GBG"]["longitude"],loca["GBG"]["latitude"]],types=gethourlyparams1())
        with open('filename.json', 'w') as f:
            json.dump(result,f)
def dumpa(name,result):
    with open(name, 'w') as f:
            json.dump(result,f)
def loada(name):
    with open(name,"r") as f:
        return json.load(f)

def checkwater(lat,lon):
    return requests.get("https://is-on-water.balbona.me/api/v1/get/"+str(lat)+"/"+str(lon)).json()
def createrandomcordinates(n=200):
    lats,lons=np.random.uniform(-90+1e-10,90,n),np.random.uniform(-180+1e-10,180,n)
    watercords,landcords=[],[]
    for i in range(n):
        print(i/n,i)
        response=checkwater(lats[i],lons[i])
        print(response)
        if response["isWater"]:
            watercords.append([lats[i],lons[i]])
        else:
            landcords.append([lats[i],lons[i]])
    dumpa("landandwatercord1000.json",{"landcords":landcords,"watercords":watercords})
def combine():
    listan=loada("landandwatercord1000.json")["landcords"]+loada("landandwatercord.json")["landcords"]
    dumpa("combineisland.json",{"cord":listan})
def createmuchdata(vikt=25):
    cordaa=loada("combineisland.json")["cord"]
    cords=cordaa[vikt:min(vikt+3,320)]
    for i in range(len(cords)):
        data=getdata(cords[i],gethourlyparams1()).json()
        dumpa(f"weatherdata//weatherdata{vikt+i}.json",{"weatherdata":data,"cords":cords[i]})
        print((vikt+i)/len(cordaa),i+vikt,len(cordaa))
    print("made")
    #dumpa("weatherdata.json",{"weatherdata":data,"cords":cords})
def createparams(data,time=[2025,4,2,0],cords=[57,12],exTime=-1):
    #print(time)
    if exTime>=14*24:
        exTime=-1
    if exTime==-1:
        time[3]=np.random.randint(0,23)
        time[2]=np.random.randint(2,15)
    else:
        time[3]=int((time[3]+exTime)%24)
        time[2]=int((time[2]+exTime//24)%28)
        #print(time[2],exTime)

    timestr=list(str(time[0])+"-04-02T00:00")
    for i in range(len(str(time[1]))):
        timestr[6-i]=str(time[1])[len(str(time[1]))-1-i]
    for i in range(len(str(time[2]))):
        timestr[9-i]=str(time[2])[len(str(time[2]))-1-i]
    for i in range(len(str(time[3]))):
        timestr[12-i]=str(time[3])[len(str(time[3]))-1-i]
    timestr="".join(timestr)
    #print(timestr,time)
    if timestr in data["hourly"]["time"]:
        indexforsearch=data["hourly"]["time"].index(timestr)
        params={"year":time[0],"month":time[1],"day":time[2],"hour":time[3],"Ta":data["hourly"]["temperature_2m"][indexforsearch],"RH":data["hourly"]["relative_humidity_2m"][indexforsearch],"Ws":data["hourly"]["wind_speed_10m"][indexforsearch],"sky":"Clear (100%)","c1":data["hourly"]["cloud_cover_low"][indexforsearch],"c2":data["hourly"]["cloud_cover_mid"][indexforsearch],"c3":data["hourly"]["cloud_cover_high"][indexforsearch], "loc":{'latitude': cords[0],'longitude': cords[1], 'altitude': data["elevation"]},"UTC":0}
        #print(params)
        return params
    #print(timestr)
    return None
def flaska(time=[2025,4,2,0]):
    varians=0
    VARIANSTMRT=0
    lll=1
    paramsfull=[None]*323*lll
    resultsreal=[None]*323*lll
    for j in range(lll):
        for i in range(323):
            data=loada(f"weatherdata//weatherdata{i}.json")
            cords=data["cords"]
            data=data["weatherdata"]
            params=createparams(data,time,cords,-1)
            paramsfull[j*323+i]=params
      #print(params)
            c1,c2,c3,t1,t2,t3=flaskapp.index(params)
            if abs(t1-t2)>1e-4:
                print(f"\t not equal {abs(t1-t2):.5f}\n{params}")
                if abs(t1-t2)>5:
                    print(f"Tmrt: \t {t1:.5f}\t{t2:.5f}\t{t3:.5f}")
                    print(i,j)
                    break
            #print(c1,c2,c3)
            print(f"Tmrt: \t {t1:.5f}\t{t2:.5f}\t{t3:.5f}")
            resultsreal[j*323+i]=c2
            varians+=(c1-c2)**2
            VARIANSTMRT+=(t1-t2)**2
        if i!=321:
            break
    print(varians,VARIANSTMRT)
    if(varians==0):
        print("\n".join(["".join([["d","a"][np.random.randint(0,2)] for i in range(np.random.randint(7,20))]) for i in range(np.random.randint(2,13))])+"\n"+"d"+"o"*np.random.randint(1,13))
    fig, ax=plt.subplots(3,3)
    dif=["Ta","RH","Ws","c1","c2","c3",'longitude', 'latitude', 'altitude']
    limess=[[-40,50],[0,100],[0,50],[0,100],[0,100],[0,100],[-180,180],[-90,90],[0,6000]]
    for ii in range(9):
        if ii>5:
            ax[ii//3][ii%3].scatter(np.array([i["loc"][dif[ii]] for i in paramsfull]),np.array(resultsreal),linewidths=0.3,alpha=0.2)
        else:
            ax[ii//3][ii%3].scatter(np.array([i[dif[ii]] for i in paramsfull]),np.array(resultsreal),linewidths=0.3,alpha=0.2)
        if ii==0:
            ax[ii//3][ii%3].plot([-40,50],[-40,50])
        ax[ii//3][ii%3].set_xlim((limess[ii][0],limess[ii][1]))
        ax[ii//3][ii%3].set_ylim((-40,60))

    plt.show()
"""
loca={"GBG":{'longitude': 12.0, 'latitude': 57.7, 'altitude': 3.}}
params={"GBG":{"year":2025,"month":1,"day":12,"hour":12,"Ta":20,"RH":50,"Ws":20,"sky":"Clear (100%)","c1":2,"c2":95,"c3":2, "loc":loca["GBG"],"UTC":1}}

flaskapp.index(params["GBG"])
"""

#createrandomcordinates(1000)
#combine()
#start=319
#createmuchdata(start)
flaska()
